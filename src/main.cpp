#include "heragrid/builder.hpp"
#include <LHAPDF/LHAPDF.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>

using namespace heragrid;
namespace fs=std::filesystem;
namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
struct PDFProvider {
    std::unique_ptr<LHAPDF::PDF> pdf;
    explicit PDFProvider(const std::string& name,int member):pdf(LHAPDF::mkPDF(name,member)) {}
    Pdf callback() {
        return [this](double x,double q) {
            if(!pdf->inRangeXQ(x,q))throw std::runtime_error("Grid node outside PDF x,Q range; extrapolation is disabled");
            PdfValues v{};
            for(int i=0;i<13;++i) {int pdg=i==6?21:i-6;if(pdf->hasFlavor(pdg))v[i]=pdf->xfxQ(pdg,x,q);}
            return v;
        };
    }
    AlphaS alpha() {return [this](double q){return pdf->alphasQ(q);};}
    std::function<double(int,double,double)> fx() {
        return [this](int slot,double x,double q) {
            if(!pdf->inRangeXQ(x,q))throw std::runtime_error("Direct convolution outside PDF x,Q range");
            int pdg=slot==6?21:slot-6;
            return pdf->hasFlavor(pdg)?pdf->xfxQ(pdg,x,q)/x:0.;
        };
    }
};
std::ofstream output(const fs::path& path) {
    if(path.has_parent_path())fs::create_directories(path.parent_path());
    std::ofstream f(path);if(!f)throw std::runtime_error("Cannot write "+path.string());f<<std::setprecision(17);return f;
}
bool validate(const Grid& grid,PDFProvider& provider,double rtol,double atol,const std::string& path,
              double reference_rtol=-1,int repeats=3) {
    if(grid.backend_digest!=backend_digest())throw std::runtime_error("Grid and direct backend source fingerprints differ; rebuild with matching physics sources");
    auto xf=provider.callback();auto alpha=provider.alpha();
    std::vector<double> fast,reference;auto start=Clock::now();
    for(int r=0;r<repeats;++r)fast=grid.convolve(xf,alpha);
    double fast_time=elapsed(start)/repeats;
    Settings reference_settings=grid.settings;
    reference_settings.integration_accuracy=reference_rtol>0?reference_rtol:std::min(grid.settings.integration_accuracy,1e-5);
    reference_settings.integration_depth=std::max(grid.settings.integration_depth,10U);
    Direct direct(reference_settings);start=Clock::now();
    auto fx=provider.fx();
    for(const auto& b:grid.bins)reference.push_back(direct.evaluate_fx(b.point,fx,alpha));
    double direct_time=elapsed(start),max_relative=0,max_scaled=0;bool ok=true;
    auto f=output(path);f<<"bin,dataset,source_row,process,effective_order,x,Q2,y,grid,direct,absolute_error,relative_error,pass\n";
    for(size_t i=0;i<fast.size();++i) {
        double err=std::abs(fast[i]-reference[i]),rel=err/std::max(std::abs(reference[i]),1e-300);
        bool pass=err<=atol+rtol*std::abs(reference[i]);ok=ok&&pass;max_relative=std::max(max_relative,rel);
        max_scaled=std::max(max_scaled,err/(atol+rtol*std::abs(reference[i])));
        const auto& b=grid.bins[i];const auto& p=b.point;
        f<<i<<','<<p.dataset<<','<<p.source_row<<','<<process_name(p.process)<<','<<b.effective_order<<','<<p.x<<','<<p.Q2<<','<<p.y
         <<','<<fast[i]<<','<<reference[i]<<','<<err<<','<<rel<<','<<pass<<'\n';
    }
    std::cout<<"Validation "<<(ok?"PASS":"FAIL")<<": max relative error="<<max_relative<<", max error/tolerance="<<max_scaled
             <<", direct="<<direct_time<<" s, grid="<<fast_time<<" s, speedup="<<direct_time/fast_time<<"\n";
    auto report=output(path+".summary.txt");report<<"pass="<<ok<<"\nmax_relative_error="<<max_relative<<"\nmax_scaled_error="<<max_scaled
        <<"\nrelative_tolerance="<<rtol<<"\nabsolute_tolerance="<<atol<<"\ndirect_seconds="<<direct_time<<"\ngrid_seconds="<<fast_time
        <<"\nspeedup="<<direct_time/fast_time<<"\npdf="<<provider.pdf->set().name()<<"\nmember="<<provider.pdf->memberID()
        <<"\nreference_integration_rtol="<<reference_settings.integration_accuracy
        <<"\nreference_integration_depth="<<reference_settings.integration_depth<<'\n';
    return ok;
}
void help() {
    std::cout<<"hera-grid generate CONFIG.yaml [--limit N]\n"
             <<"hera-grid inspect GRID.hgrid\n"
             <<"hera-grid evaluate GRID.hgrid --pdf SET [--member N] [--output CSV] [--order 1|2]\n"
             <<"hera-grid validate GRID.hgrid --pdf SET [--member N] [--output CSV] [--rtol R] [--atol A] [--reference-rtol R]\n"
             <<"Generation validates against direct CTEQDisPac before saving; failure exits 2.\n";
}
}
int main(int argc,char** argv) {
    try {
        if(argc<2 || std::string(argv[1])=="--help") {help();return 0;}
        if(argc<3)throw std::runtime_error("Missing configuration/grid path");
        std::string command=argv[1],path=argv[2],pdf,output_path;int member=0,order=-1;size_t limit=0;
        std::set<std::string> allowed;
        if(command=="generate")allowed={"--limit"};
        else if(command=="evaluate")allowed={"--pdf","--member","--output","--order"};
        else if(command=="validate")allowed={"--pdf","--member","--output","--rtol","--atol","--reference-rtol"};
        else if(command!="inspect")throw std::runtime_error("Unknown command "+command);
        double rtol=.001,atol=1e-8,reference_rtol=-1;
        for(int i=3;i<argc;++i) {
            std::string flag=argv[i];
            if(!allowed.count(flag))throw std::runtime_error("Option "+flag+" is not valid for "+command);
            if(i+1>=argc)throw std::runtime_error("Missing value for "+flag);
            std::string value=argv[++i];
            if(flag=="--pdf")pdf=value;else if(flag=="--member")member=std::stoi(value);
            else if(flag=="--output")output_path=value;else if(flag=="--order")order=std::stoi(value);
            else if(flag=="--limit") {
                if(value.empty() || value[0]=='-')throw std::runtime_error("Invalid --limit");
                limit=std::stoull(value);
                if(!limit)throw std::runtime_error("--limit must be positive");
            }else if(flag=="--rtol")rtol=std::stod(value);
            else if(flag=="--atol")atol=std::stod(value);
            else if(flag=="--reference-rtol")reference_rtol=std::stod(value);
            else throw std::runtime_error("Unknown option "+flag);
        }
        if(member<0 || !(rtol>0) || !(atol>=0) || !std::isfinite(rtol) || !std::isfinite(atol))throw std::runtime_error("Invalid arguments");
        if(command=="generate") {
            Config config=read_config(path);
            if(limit && limit<config.points.size())config.points.resize(limit);
            for(int o:config.orders) {
                std::string dest=config.output+(o==1?"_NLO.hgrid":"_NNLO.hgrid");
                if(fs::exists(dest))throw std::runtime_error("Output already exists: "+dest);
            }
            std::cout<<"Generating PDF-independent coefficient grid for "<<config.points.size()<<" points; "<<config.settings.nodes<<" x nodes\n";
            if(config.settings.order==2 && config.settings.allow_nlo_cc)
                std::cout<<"Source convention: NC includes NNLO photon terms; CC remains NLO (stored per bin).\n";
            PDFProvider provider(config.pdf,config.member); // fail early if validation PDF unavailable
            auto start=Clock::now();
            Grid full=build(config.points,config.settings,[&](size_t done,size_t total) {
                std::cout<<"Completed "<<done<<"/"<<total<<" ("<<elapsed(start)<<" s)"<<std::endl;
            },config.work_directory);
            bool all_ok=true;
            for(int o:config.orders) {
                Grid grid=full;grid.settings.order=o;
                for(auto& b:grid.bins) {b.weights.resize((o+1)*13*grid.settings.nodes);b.effective_order=std::min(o,b.effective_order);}
                std::string dest=config.output+(o==1?"_NLO.hgrid":"_NNLO.hgrid");
                std::string report=(fs::path(config.work_directory)/(fs::path(dest).filename().string()+".validation.csv")).string();
                bool ok=validate(grid,provider,config.validation_relative,config.validation_absolute,report);
                if(ok) {grid.save(dest);std::cout<<"Saved "<<dest<<'\n';}
                else {all_ok=false;std::cerr<<"Grid not saved: refine XNodes/integration or examine validation report\n";}
            }
            return all_ok?0:2;
        }
        Grid grid=Grid::load(path);
        if(command=="inspect") {
            std::cout<<"Format: HERAGRID 1\nBackend: CTEQDisPac (SACOT-chi source implementation)\nSource SHA256: "<<grid.backend_digest
                     <<"\nRequested order: "<<(grid.settings.order==1?"NLO":"NNLO")<<"\nBins: "<<grid.bins.size()<<"\nNodes/bin: "<<grid.settings.nodes<<'\n';
            size_t cc=0;for(auto& b:grid.bins)if(b.point.process==Process::CCem || b.point.process==Process::CCep)++cc;
            std::cout<<"CC bins (NLO only): "<<cc<<"\nScale: muF=muR=Q, exact per-bin Q\n";return 0;
        }
        if(pdf.empty())throw std::runtime_error("--pdf is required");
        PDFProvider provider(pdf,member);
        if(command=="validate") {
            if(order!=-1 || limit)throw std::runtime_error("validate uses all bins at stored order");
            if(reference_rtol!=-1 && !(reference_rtol>0 && reference_rtol<=.01))throw std::runtime_error("Invalid --reference-rtol");
            return validate(grid,provider,rtol,atol,output_path.empty()?path+".check.csv":output_path,reference_rtol)?0:2;
        }
        if(command=="evaluate") {
            auto values=grid.convolve(provider.callback(),provider.alpha(),order);
            auto f=output(output_path.empty()?path+".predictions.csv":output_path);
            f<<"bin,dataset,source_row,process,x,Q2,y,prediction\n";
            for(size_t i=0;i<values.size();++i) {
                const auto& p=grid.bins[i].point;
                f<<i<<','<<p.dataset<<','<<p.source_row<<','<<process_name(p.process)<<','<<p.x<<','<<p.Q2<<','<<p.y<<','<<values[i]<<'\n';
            }
            std::cout<<"Evaluated "<<values.size()<<" bins with "<<pdf<<"/"<<member<<'\n';return 0;
        }
        throw std::runtime_error("Unknown command "+command);
    }catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
