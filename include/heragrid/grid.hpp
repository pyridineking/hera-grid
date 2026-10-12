#pragma once
#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace heragrid {
// LHAPDF convention for VALUES: x*f. Slot = signed PDG + 6, gluon in slot 6.
using PdfValues = std::array<double,13>;
using Pdf = std::function<PdfValues(double x,double Q)>;
using AlphaS = std::function<double(double Q)>;
enum class Process { NCep=0, NCem=1, CCep=2, CCem=3 };
std::string process_name(Process p);
Process process_from_name(const std::string& s);
struct Point {
    std::string dataset;
    std::string source;
    std::size_t source_row=0; // 0-based data row, excluding headers
    Process process=Process::NCep;
    int flavor=0; // 0=inclusive, 14=charm, 15=bottom
    bool reduced=true;
    double x=0,Q2=0,y=0;
};
void validate_point(const Point& p);
struct Settings {
    int order=2; // maximum alpha_s power, 1=NLO, 2=NNLO
    int nodes=64;
    int threads=1;
    double mc=1.3,mb=4.75,mt=1e99;
    double mw=80.835,gf=1.16638e-5;
    double integration_accuracy=1e-4;
    unsigned integration_depth=8;
    // A flattening factor for xf, independent of any reference PDF.
    double small_x_power=-.3,large_x_power=1.;
    // Source backend has only NLO CC. Must explicitly acknowledge this for NNLO+CC.
    bool allow_nlo_cc=false;
};
void validate_settings(const Settings& s);
struct Bin {
    Point point;
    int effective_order=0;
    std::vector<double> nodes;
    // Layout [power][parton slot][x node]. Powers multiply alpha_s(Q)^power, NOT alpha_s/(2*pi).
    std::vector<double> weights;
};
class Grid {
public:
    Settings settings;
    std::string backend_digest;
    std::vector<Bin> bins;
    void save(const std::string& path) const;
    static Grid load(const std::string& path);
    // Each row is at exactly the experimental x,Q2,y. Interpolation is in convolution x.
    std::vector<double> convolve(const Pdf& xf,const AlphaS& alpha,int order=-1) const;
    double convolve_bin(std::size_t bin,const Pdf& xf,const AlphaS& alpha,int order=-1) const;
    void check() const;
};
double envelope(double x,const Settings& s);
std::vector<double> chebyshev_nodes(double xmin,int n);
double chebyshev_basis(double x,double xmin,int mode,const Settings& s);
// Source-compatible reduced HERA cross-section and CC non-reduced normalization.
double cross_section(const Point& p,const std::vector<double>& sf,const Settings& s);
} // namespace heragrid
