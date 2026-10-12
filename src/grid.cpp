#include "heragrid/grid.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace heragrid {
namespace {
constexpr double pi=3.14159265358979323846;
void require(bool b,const std::string& msg) {if(!b)throw std::runtime_error(msg);}
template<class T> void read(std::istream& in,T& v) {require(bool(in>>v),"Truncated or invalid HERA grid");}
}
std::string process_name(Process p) {
    switch(p) {case Process::NCep:return "NCep";case Process::NCem:return "NCem";
               case Process::CCep:return "CCep";case Process::CCem:return "CCem";}
    throw std::runtime_error("Invalid process");
}
Process process_from_name(const std::string& s) {
    if(s=="NCep")return Process::NCep;
    if(s=="NCem")return Process::NCem;
    if(s=="CCep")return Process::CCep;
    if(s=="CCem")return Process::CCem;
    throw std::runtime_error("Unsupported HERA process: "+s);
}
void validate_point(const Point& p) {
    process_name(p.process);
    require(p.x>0 && p.x<1 && p.Q2>0 && std::isfinite(p.Q2) && p.y>=0 && p.y<=1,"Invalid HERA x,Q2,y");
    require(p.flavor==0 || p.flavor==14 || p.flavor==15,"Flavor must be 0,14,15");
    bool nc=p.process==Process::NCep || p.process==Process::NCem;
    require(nc || p.flavor==0,"Heavy-tagged CC not supported");
    require(!nc || p.reduced,"Non-reduced NC is rejected: original Hera.cxx applies a CC prefactor here");
}
void validate_settings(const Settings& s) {
    require(s.order==1 || s.order==2,"Order must be NLO (1) or NNLO (2)");
    require(s.nodes>=8 && s.nodes<=512,"XNodes must be 8..512");
    require(s.threads>=1 && s.threads<=256,"Threads must be 1..256");
    require(s.mc>.2 && s.mb>s.mc && s.mt>s.mb && std::isfinite(s.mt),"Heavy-quark masses must lie above source strange mass 0.2 and be ordered");
    require(s.mw>0 && s.gf>0 && std::isfinite(s.mw) && std::isfinite(s.gf),"Invalid EW constants");
    require(s.integration_accuracy>0 && s.integration_accuracy<=.01 && s.integration_depth<=25,"Invalid integration settings");
    require(std::isfinite(s.small_x_power) && std::abs(s.small_x_power)<=2 &&
            s.large_x_power>=0 && s.large_x_power<=5,"Invalid PDF interpolation envelope");
}
double envelope(double x,const Settings& s) {
    return std::pow(x,s.small_x_power)*std::pow(1-x,s.large_x_power);
}
std::vector<double> chebyshev_nodes(double xmin,int n) {
    std::vector<double> nodes(n);
    for(int j=0;j<n;++j) {
        double t=std::cos(pi*(j+.5)/n);
        nodes[j]=std::exp(.5*(1-t)*std::log(xmin));
    }
    return nodes; // descending x, no x=1 endpoint
}
double chebyshev_basis(double x,double xmin,int mode,const Settings& s) {
    if(x>=1)return 0;
    require(x>=xmin*(1-1e-10) && x>0,"Backend requested a PDF below this bin's convolution domain");
    double t=1-2*std::log(x)/std::log(xmin);
    t=std::clamp(t,-1.,1.);
    return envelope(x,s)*std::cos(mode*std::acos(t));
}
double cross_section(const Point& p,const std::vector<double>& f,const Settings& s) {
    double yp=1+std::pow(1-p.y,2),ym=1-std::pow(1-p.y,2),v=0;
    switch(p.process) {
        case Process::NCep:v=f.at(0)-ym/yp*f.at(2)-p.y*p.y/yp*f.at(1);break;
        case Process::NCem:v=f.at(0)+ym/yp*f.at(2)-p.y*p.y/yp*f.at(1);break;
        case Process::CCep:v=.5*(yp*f.at(0)-ym*f.at(2)-p.y*p.y*f.at(1));break;
        case Process::CCem:v=.5*(yp*f.at(0)+ym*f.at(2)-p.y*p.y*f.at(1));break;
    }
    if(!p.reduced)v*=s.gf*s.gf/(2*3.1415927*p.x)*std::pow(1+p.Q2/(s.mw*s.mw),-2)*.389379338e9;
    require(std::isfinite(v),"Non-finite structure function/cross-section");return v;
}
void Grid::check() const {
    validate_settings(settings);
    require(!bins.empty() && bins.size()<=100000,"Invalid grid bin count");
    for(const auto& b:bins) {
        validate_point(b.point);
        bool cc=b.point.process==Process::CCep || b.point.process==Process::CCem;
        require(b.effective_order==(cc?1:settings.order),"Inconsistent perturbative order metadata");
        require(!(cc && settings.order==2) || settings.allow_nlo_cc,"NNLO grid contains unacknowledged NLO-only CC");
        require(b.nodes.size()==size_t(settings.nodes),"Invalid x node count");
        require(b.weights.size()==size_t((settings.order+1)*13*settings.nodes),"Invalid weight count");
        for(size_t j=0;j<b.nodes.size();++j) {
            require(std::isfinite(b.nodes[j]) && b.nodes[j]>b.point.x && b.nodes[j]<1,"Invalid convolution x node");
            if(j)require(b.nodes[j]<b.nodes[j-1],"Grid nodes not strictly descending");
        }
        for(double w:b.weights)require(std::isfinite(w),"Non-finite grid weight");
        if(cc && settings.order==2)for(size_t j=26*b.nodes.size();j<b.weights.size();++j)
            require(b.weights[j]==0,"CC NNLO coefficients must be zero for this backend");
    }
}
void Grid::save(const std::string& path) const {
    check();
    std::filesystem::path dest(path);
    if(dest.has_parent_path())std::filesystem::create_directories(dest.parent_path());
    // Refuse overwrites so an interrupted generation cannot corrupt an existing grid.
    require(!std::filesystem::exists(dest),"Output already exists; choose a new grid filename: "+path);
    std::ofstream out(path);require(bool(out),"Cannot write "+path);out<<std::setprecision(17);
    const auto& s=settings;
    out<<"HERAGRID 1\n"<<std::quoted(backend_digest)<<"\n";
    out<<s.order<<' '<<s.nodes<<' '<<s.mc<<' '<<s.mb<<' '<<s.mt<<' '<<s.mw<<' '<<s.gf<<' '
       <<s.integration_accuracy<<' '<<s.integration_depth<<' '<<s.small_x_power<<' '<<s.large_x_power<<' '<<s.allow_nlo_cc<<' '<<bins.size()<<'\n';
    for(const auto& b:bins) {
        const auto& p=b.point;
        out<<std::quoted(p.dataset)<<' '<<std::quoted(p.source)<<' '<<p.source_row<<' '<<process_name(p.process)<<' '
           <<p.flavor<<' '<<p.reduced<<' '<<p.x<<' '<<p.Q2<<' '<<p.y<<' '<<b.effective_order<<'\n';
        for(double x:b.nodes)out<<x<<' ';
        out<<'\n';
        for(double w:b.weights)out<<w<<' ';
        out<<'\n';
    }
    out<<"END_HERAGRID\n";out.flush();require(bool(out),"Grid write failed: "+path);
}
Grid Grid::load(const std::string& path) {
    std::ifstream in(path);require(bool(in),"Cannot read grid "+path);
    std::string magic;int version=0;read(in,magic);read(in,version);
    require(magic=="HERAGRID" && version==1,"Unsupported HERA grid format");
    Grid g;auto& s=g.settings;require(bool(in>>std::quoted(g.backend_digest)),"Missing backend provenance");
    read(in,s.order);read(in,s.nodes);read(in,s.mc);read(in,s.mb);read(in,s.mt);read(in,s.mw);read(in,s.gf);
    read(in,s.integration_accuracy);read(in,s.integration_depth);read(in,s.small_x_power);read(in,s.large_x_power);read(in,s.allow_nlo_cc);
    validate_settings(s);size_t count=0;read(in,count);require(count>0 && count<=100000,"Invalid bin count");
    for(size_t i=0;i<count;++i) {
        Bin b;auto& p=b.point;std::string process;
        require(bool(in>>std::quoted(p.dataset)>>std::quoted(p.source)),"Invalid bin labels");
        read(in,p.source_row);read(in,process);p.process=process_from_name(process);
        read(in,p.flavor);read(in,p.reduced);read(in,p.x);read(in,p.Q2);read(in,p.y);read(in,b.effective_order);
        b.nodes.resize(s.nodes);b.weights.resize((s.order+1)*13*s.nodes);
        for(double& x:b.nodes)read(in,x);
        for(double& w:b.weights)read(in,w);
        g.bins.push_back(std::move(b));
    }
    read(in,magic);require(magic=="END_HERAGRID","Incomplete grid (missing footer)");
    require(!(in>>magic),"Unexpected data after grid footer");g.check();return g;
}
std::vector<double> Grid::convolve(const Pdf& xf,const AlphaS& alpha,int order) const {
    std::vector<double> result;result.reserve(bins.size());
    for(size_t i=0;i<bins.size();++i)result.push_back(convolve_bin(i,xf,alpha,order));
    return result;
}
double Grid::convolve_bin(size_t bin,const Pdf& xf,const AlphaS& alpha,int order) const {
    if(order<0)order=settings.order;
    require(order>=0 && order<=settings.order,"Requested order exceeds grid order");
    const auto& b=bins.at(bin);
        const size_t n=b.nodes.size();double q=std::sqrt(b.point.Q2),as=alpha(q);
        require(std::isfinite(as) && as>=0,"Invalid alpha_s callback");
        long double sum=0;
        for(size_t j=0;j<n;++j) {
            auto values=xf(b.nodes[j],q);
            for(int f=0;f<13;++f) {
                require(std::isfinite(values[f]),"Non-finite xf callback");
                double power=1;
                for(int p=0;p<=order;++p) {
                    sum+=static_cast<long double>(power)*b.weights[(p*13+f)*n+j]*values[f];power*=as;
                }
            }
        }
        require(std::isfinite(double(sum)),"Non-finite convolution result");
    return double(sum);
}
} // namespace heragrid
