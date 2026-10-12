#include "heragrid/builder.hpp"
#include "CTEQDisPac.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace heragrid {
std::string backend_digest() { return HERA_BACKEND_DIGEST; }
namespace {
int interaction(Process p) {
    if(p==Process::CCep)return CTEQDisPac::kWp;
    if(p==Process::CCem)return CTEQDisPac::kWm;
    return CTEQDisPac::kCombinedNC;
}
Bin build_bin(const Point& point,const Settings& settings) {
    validate_point(point);
    const bool cc=point.process==Process::CCep || point.process==Process::CCem;
    if(cc && settings.order==2 && !settings.allow_nlo_cc)
        throw std::runtime_error("CTEQDisPac CC is NLO-only. Set AllowNLOChargedCurrent: true to explicitly retain source behavior.");
    Bin bin;bin.point=point;bin.effective_order=cc?1:settings.order;
    bin.nodes=chebyshev_nodes(point.x,settings.nodes);
    const int n=settings.nodes;bin.weights.assign((settings.order+1)*13*n,0);
    CTEQDisPac backend;
    backend.SetHFMass(settings.mc,settings.mb,settings.mt);
    backend.SetGridIntegration(settings.integration_accuracy,settings.integration_depth);
    const double q=std::sqrt(point.Q2),probe_alpha=.2;
    std::vector<double> modes((settings.order+1)*n);
    // Global Chebyshev modes avoid compact-support basis functions being missed by adaptive quadrature.
    for(int f=0;f<13;++f) {
        // The backend never samples flavors above its factorization-scale threshold.
        int flavor=std::abs(f-6);
        if((flavor==4 && q<=settings.mc) || (flavor==5 && q<=settings.mb) || (flavor==6 && q<=settings.mt))continue;
        for(int mode=0;mode<n;++mode) {
            auto pdf=[&](int slot,double x,double mu) {
                if(slot!=f)return 0.;
                if(std::abs(mu/q-1)>1e-10)throw std::runtime_error("Backend PDF scale is not Q; regenerate with a scale-aware format");
                return chebyshev_basis(x,point.x,mode,settings)/x;
            };
            auto eval=[&](int order,double alpha) {
                backend.SetOrder(order);
                backend.SetPDFAlphaS(pdf,[=](double){return alpha;});
                return cross_section(point,backend.CalculateSF(interaction(point.process),point.flavor,point.x,q),settings);
            };
            double lo=eval(1,0),nlo=eval(1,probe_alpha);
            modes[mode]=lo;modes[n+mode]=(nlo-lo)/probe_alpha;
            if(settings.order==2)modes[2*n+mode]=cc?0:(eval(2,probe_alpha)-nlo)/(probe_alpha*probe_alpha);
        }
        // DCT transforms responses to polynomial modes into reusable weights on x*f at nodes.
        for(int p=0;p<=settings.order;++p)for(int j=0;j<n;++j) {
            long double w=modes[p*n];
            for(int mode=1;mode<n;++mode)
                w+=2*modes[p*n+mode]*std::cos(mode*std::acos(-1.)*(j+.5)/n);
            bin.weights[(p*13+f)*n+j]=double(w/n)/envelope(bin.nodes[j],settings);
        }
    }
    return bin;
}
}
Direct::Direct(const Settings& s):settings_(s),backend_(std::make_unique<CTEQDisPac>()) {
    validate_settings(s);backend_->SetHFMass(s.mc,s.mb,s.mt);
    backend_->SetGridIntegration(s.integration_accuracy,s.integration_depth);
}
Direct::~Direct()=default;
double Direct::evaluate(const Point& p,const Pdf& xf,const AlphaS& alpha,int order) {
    return evaluate_fx(p,[&](int f,double x,double q){return xf(x,q).at(f)/x;},alpha,order);
}
double Direct::evaluate_fx(const Point& p,const std::function<double(int,double,double)>& fx,const AlphaS& alpha,int order) {
    validate_point(p);if(order<0)order=settings_.order;
    if(order!=1 && order!=2)throw std::runtime_error("Direct order must be NLO or NNLO");
    backend_->SetOrder(order);
    backend_->SetPDFAlphaS([&](int f,double x,double q){return x>=1?0:fx(f,x,q);},alpha);
    return cross_section(p,backend_->CalculateSF(interaction(p.process),p.flavor,p.x,std::sqrt(p.Q2)),settings_);
}
Grid build(const std::vector<Point>& points,const Settings& s,
           const std::function<void(size_t,size_t)>& progress,const std::string& checkpoint_directory) {
    validate_settings(s);if(points.empty())throw std::runtime_error("No points selected");
    for(const auto& p:points) {
        validate_point(p);
        if(s.order==2 && !s.allow_nlo_cc && (p.process==Process::CCep || p.process==Process::CCem))
            throw std::runtime_error("NNLO request contains CC: explicitly set AllowNLOChargedCurrent: true (backend has only NLO CC)");
    }
    Grid grid;grid.settings=s;grid.backend_digest=backend_digest();grid.bins.resize(points.size());
    namespace fs=std::filesystem;
    if(!checkpoint_directory.empty())fs::create_directories(checkpoint_directory);
    std::atomic<size_t> next{0};std::atomic<bool> failed{false};std::mutex mutex;std::exception_ptr error;size_t done=0;
    auto worker=[&] {
        try {
            while(!failed) {
                size_t i=next++;if(i>=points.size())break;
                fs::path checkpoint=fs::path(checkpoint_directory)/("bin_"+std::to_string(i)+".hgrid");
                if(!checkpoint_directory.empty() && fs::exists(checkpoint)) {
                    Grid saved=Grid::load(checkpoint.string());const auto& t=saved.settings;const auto& p=points[i];
                    if(saved.bins.size()!=1)throw std::runtime_error("Invalid checkpoint bin count");
                    const auto& q=saved.bins[0].point;
                    if(saved.backend_digest!=grid.backend_digest || t.order!=s.order || t.nodes!=s.nodes ||
                       t.mc!=s.mc || t.mb!=s.mb || t.mt!=s.mt || t.mw!=s.mw || t.gf!=s.gf ||
                       t.integration_accuracy!=s.integration_accuracy || t.integration_depth!=s.integration_depth ||
                       t.small_x_power!=s.small_x_power || t.large_x_power!=s.large_x_power || t.allow_nlo_cc!=s.allow_nlo_cc ||
                       p.dataset!=q.dataset || p.source!=q.source || p.source_row!=q.source_row ||
                       p.process!=q.process || p.flavor!=q.flavor || p.reduced!=q.reduced || p.x!=q.x || p.Q2!=q.Q2 || p.y!=q.y)
                        throw std::runtime_error("Checkpoint settings/point mismatch; choose a new Output prefix");
                    grid.bins[i]=std::move(saved.bins[0]);
                } else {
                    grid.bins[i]=build_bin(points[i],s);
                    if(!checkpoint_directory.empty()) {
                        Grid saved;saved.settings=s;saved.backend_digest=grid.backend_digest;saved.bins={grid.bins[i]};
                        fs::path partial=checkpoint.string()+".partial";
                        // This is exclusively our known incomplete checkpoint, never an arbitrary runtime file.
                        if(fs::exists(partial))fs::remove(partial);
                        saved.save(partial.string());fs::rename(partial,checkpoint);
                    }
                }
                std::lock_guard<std::mutex> lock(mutex);++done;if(progress)progress(done,points.size());
            }
        }catch(...) {
            std::lock_guard<std::mutex> lock(mutex);if(!error)error=std::current_exception();failed=true;
        }
    };
    std::vector<std::thread> threads;
    for(int i=0;i<std::min<int>(s.threads,points.size());++i)threads.emplace_back(worker);
    for(auto& t:threads)t.join();
    if(error)std::rethrow_exception(error);
    grid.check();return grid;
}
} // namespace heragrid
