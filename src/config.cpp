#include "heragrid/builder.hpp"
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace heragrid {
namespace fs=std::filesystem;
namespace {
void require(bool b,const std::string& s) {if(!b)throw std::runtime_error(s);}
std::string trim(std::string s) {
    auto a=s.find_first_not_of(" \r\n\t"),b=s.find_last_not_of(" \r\n\t");
    return a==std::string::npos ? "" : s.substr(a,b-a+1);
}
fs::path resolve(const fs::path& base,const std::string& p) {
    fs::path out(p);return fs::absolute(out.is_absolute()?out:base/out).lexically_normal();
}
int parse_order(const YAML::Node& n) {
    std::string s=n.as<std::string>();if(s=="NLO" || s=="1")return 1;if(s=="NNLO" || s=="2")return 2;
    throw std::runtime_error("Order must be NLO or NNLO, received "+s);
}
void load_dataset(const fs::path& cfg,const fs::path& data_root,std::vector<Point>& points) {
    YAML::Node ds=YAML::LoadFile(cfg.string());
    const std::string cls=ds["Class"].as<std::string>();
    require(cls=="HeraData" || cls=="HeraDataNLO","Only HeraData/HeraDataNLO dataset configs supported: "+cfg.string());
    std::string theory=ds["Theory"].as<std::string>();
    require(theory=="Internal" || theory=="Internal_c" || theory=="Internal_b","Unsupported HERA theory "+theory);
    int flavor=theory=="Internal_c"?14:theory=="Internal_b"?15:0;
    std::string name=ds["Name"].as<std::string>(),sub=ds["SubDir"].as<std::string>();
    for(const auto& file:ds["DataFiles"]) {
        fs::path path=data_root/sub/file["Filename"].as<std::string>();
        std::ifstream in(path);require(bool(in),"Cannot open HERA data: "+path.string());
        auto meta=file["Metadata"];
        Process process=process_from_name(meta["Type"].as<std::string>());
        bool reduced=meta["reduced"].as<int>()!=0;
        std::vector<std::string> names;int declared_rows=-1,columns=-1;size_t row=0;
        std::string line;
        while(std::getline(in,line)) {
            line=trim(line.substr(0,line.find('#')));if(line.empty())continue;
            auto equals=line.find('=');
            if(equals!=std::string::npos) {
                std::string key=trim(line.substr(0,equals)),value=trim(line.substr(equals+1));
                if(key=="NData")declared_rows=std::stoi(value);
                else if(key=="NColumn")columns=std::stoi(value);
                else if(key=="ColumnNames") {
                    std::istringstream list(value);std::string token;
                    while(std::getline(list,token,','))names.push_back(trim(token));
                }else throw std::runtime_error("Unknown data header "+key);
                continue;
            }
            require(!names.empty() && columns==int(names.size()),"Invalid HERA column header in "+path.string());
            std::istringstream data(line);std::vector<double> values;double value;
            while(data>>value) {require(std::isfinite(value),"Non-finite input data");values.push_back(value);}
            require(data.eof() && values.size()==names.size(),"Column count mismatch in "+path.string()+" row "+std::to_string(row));
            auto get=[&](const std::string& key) {
                auto i=std::find(names.begin(),names.end(),key);require(i!=names.end(),"Missing data column "+key);
                return values[size_t(i-names.begin())];
            };
            bool selected=true;
            for(const auto& cut:ds["Cuts"]) {
                require(!cut["Formula"] || cut["Formula"].as<std::string>().empty(),"Formula cuts must be evaluated explicitly before grid generation");
                double v=get(cut["VarName"].as<std::string>());
                if(v<cut["Min"].as<double>() || v>cut["Max"].as<double>())selected=false;
            }
            if(selected) {
                Point p{name,path.string(),row,process,flavor,reduced,get("x"),get("Q2"),get("y")};
                validate_point(p);points.push_back(p);
            }
            ++row;
        }
        require(declared_rows>=0 && row==size_t(declared_rows),"NData mismatch in "+path.string());
    }
}
}
Config read_config(const std::string& path) {
    YAML::Node cfg=YAML::LoadFile(path);Config c;auto& s=c.settings;
    const std::set<std::string> allowed{"Orders","Order","Output","WorkDirectory","XNodes","Threads","Masses","Integration","Envelope",
        "AllowNLOChargedCurrent","Points","Datasets","DataPath","Validation","Electroweak","Backend"};
    require(cfg.IsMap(),"Config must be a YAML mapping");
    for(auto entry:cfg)require(allowed.count(entry.first.as<std::string>()),"Unknown config field "+entry.first.as<std::string>());
    if(cfg["Backend"])require(cfg["Backend"].as<std::string>()=="CTEQDisPac","This grid builder adapts the CTEQDisPac backend");
    require(!(cfg["Order"] && cfg["Orders"]),"Specify Order or Orders, not both");
    if(cfg["Orders"])for(auto o:cfg["Orders"])c.orders.push_back(parse_order(o));
    else if(cfg["Order"])c.orders.push_back(parse_order(cfg["Order"]));
    else throw std::runtime_error("Config must explicitly specify Order or Orders");
    std::sort(c.orders.begin(),c.orders.end());c.orders.erase(std::unique(c.orders.begin(),c.orders.end()),c.orders.end());
    require(!c.orders.empty(),"Orders list is empty");s.order=c.orders.back();
    if(cfg["XNodes"])s.nodes=cfg["XNodes"].as<int>();
    if(cfg["Threads"])s.threads=cfg["Threads"].as<int>();
    if(cfg["AllowNLOChargedCurrent"])s.allow_nlo_cc=cfg["AllowNLOChargedCurrent"].as<bool>();
    if(auto m=cfg["Masses"]) {if(m["Charm"])s.mc=m["Charm"].as<double>();if(m["Bottom"])s.mb=m["Bottom"].as<double>();if(m["Top"])s.mt=m["Top"].as<double>();}
    if(auto m=cfg["Electroweak"]) {if(m["MW"])s.mw=m["MW"].as<double>();if(m["GF"])s.gf=m["GF"].as<double>();}
    if(auto v=cfg["Integration"]) {if(v["RelativeTolerance"])s.integration_accuracy=v["RelativeTolerance"].as<double>();if(v["MaxDepth"])s.integration_depth=v["MaxDepth"].as<unsigned>();}
    if(auto v=cfg["Envelope"]) {if(v["SmallXPower"])s.small_x_power=v["SmallXPower"].as<double>();if(v["LargeXPower"])s.large_x_power=v["LargeXPower"].as<double>();}
    if(auto v=cfg["Validation"]) {
        if(v["PDF"])c.pdf=v["PDF"].as<std::string>();
        if(v["Member"])c.member=v["Member"].as<int>();
        if(v["RelativeTolerance"])c.validation_relative=v["RelativeTolerance"].as<double>();
        if(v["AbsoluteTolerance"])c.validation_absolute=v["AbsoluteTolerance"].as<double>();
    }
    require(c.member>=0 && c.validation_relative>0 && c.validation_absolute>=0 && std::isfinite(c.validation_relative) && std::isfinite(c.validation_absolute),"Invalid validation settings");
    fs::path base=fs::absolute(path).parent_path();
    c.output=resolve(base,cfg["Output"].as<std::string>()).string();
    c.work_directory=cfg["WorkDirectory"]?resolve(base,cfg["WorkDirectory"].as<std::string>()).string():c.output+".work";
    if(cfg["Datasets"]) {
        require(bool(cfg["DataPath"]),"Datasets require DataPath");
        fs::path data=resolve(base,cfg["DataPath"].as<std::string>());
        for(auto d:cfg["Datasets"])load_dataset(resolve(base,d.as<std::string>()),data,c.points);
    }
    for(const auto& n:cfg["Points"]) {
        Point p;p.dataset=n["Dataset"]?n["Dataset"].as<std::string>():"inline";
        p.source=fs::absolute(path).string();p.source_row=c.points.size();
        p.process=process_from_name(n["Process"].as<std::string>());
        p.flavor=n["Flavor"]?n["Flavor"].as<int>():0;p.reduced=n["Reduced"]?n["Reduced"].as<bool>():true;
        p.x=n["x"].as<double>();p.Q2=n["Q2"].as<double>();p.y=n["y"].as<double>();validate_point(p);c.points.push_back(p);
    }
    validate_settings(s);require(!c.points.empty(),"Config selects zero data points");return c;
}
} // namespace heragrid
