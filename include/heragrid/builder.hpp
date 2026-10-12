#pragma once
#include "heragrid/grid.hpp"
#include <memory>
class CTEQDisPac;
namespace heragrid {
std::string backend_digest();
class Direct {
public:
    explicit Direct(const Settings& settings);
    ~Direct();
    double evaluate(const Point& point,const Pdf& xf,const AlphaS& alpha,int order=-1);
    // Original fitting-package convention: f (not xf), slot - 6 is signed PDG, 6 is gluon.
    double evaluate_fx(const Point& point,const std::function<double(int,double,double)>& fx,
                       const AlphaS& alpha,int order=-1);
private:
    Settings settings_;
    std::unique_ptr<CTEQDisPac> backend_;
};
Grid build(const std::vector<Point>& points,const Settings& settings,
           const std::function<void(std::size_t,std::size_t)>& progress={},
           const std::string& checkpoint_directory="");
struct Config {
    Settings settings;
    std::vector<Point> points;
    std::vector<int> orders;
    std::string output;
    std::string work_directory;
    std::string pdf="CT18NNLO";
    int member=0;
    double validation_relative=1e-3,validation_absolute=1e-8;
};
Config read_config(const std::string& path);
} // namespace heragrid
