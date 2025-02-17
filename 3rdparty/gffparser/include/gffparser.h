#ifndef GFFPARSER_H
#define GFFPARSER_H
#include <deque>
#include <mutex>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>
#include <thread>

namespace gffparser {

enum class gff_field_type_t: int {
    SEQID = 0, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES, FIELDSLEN
};

struct gff_attr_t {
    bool is_string() const;
    bool is_integer() const;
    bool is_float() const;
    gff_attr_t& set_string(const std::string &val);
    gff_attr_t& set_integer(int64_t val);
    gff_attr_t& set_float(double val);
    const std::string& get_string() const noexcept(false);
    const int64_t& get_integer() const noexcept(false);
    const double& get_float() const noexcept(false);

    bool empty() const;
    std::string str() const;
    bool eq(const gff_attr_t &v) const;

    gff_attr_t() {}
    gff_attr_t(const std::string &val): _value(val) {}
    gff_attr_t(int64_t val): _value(val) {}
    gff_attr_t(double val): _value(val) {}
protected:
    std::variant<std::string, int64_t, double> _value;
};

struct gff_attribute_t
{
    std::string name;
    gff_attr_t value;
    std::vector<gff_attr_t> orvalues;
    gff_attribute_t(const std::string &nm): name(nm) {}
    gff_attribute_t(const std::string &nm, const std::string &val): name(nm), value(val) {}
    gff_attribute_t(const std::string &nm, int64_t val): name(nm), value(val) {}
    gff_attribute_t(const std::string &nm, double val): name(nm), value(val) {}
    gff_attribute_t& add_value(const std::string &val);
    gff_attribute_t& add_value(int64_t val);
    gff_attribute_t& add_value(double val);
};

struct gff_postition_t {
    std::string seqid;
    uint64_t start;
    uint64_t end;
    gff_postition_t(const std::string &sid, uint64_t spos, uint64_t epos)
        : seqid(sid), start(spos), end(epos) {}
    gff_postition_t(const std::string &sid, uint64_t pos)
        : seqid(sid), start(pos), end(pos) {}
    gff_postition_t(): gff_postition_t("", 0, 0) {}
    bool empty() const;
    bool singlepos() const;
    bool operator==(const gff_postition_t &pos) const;
    gff_postition_t& operator=(const gff_postition_t &pos);
    bool in(uint64_t pos) const;
    bool in(const std::string &sid, uint64_t pos) const;
    bool intersect(uint64_t spos, uint64_t epos) const;
    bool intersect(const std::string &sid, uint64_t spos, uint64_t epos) const;
};

struct gff_data_t {
    static double d_nodata;
    std::string source;
    std::string type;
    gff_postition_t position;
    double score = d_nodata;
    char strand = 0;
    char phase = 0;
    uint64_t linenum;
    std::unordered_map<std::string, gff_attr_t> attributes;
    bool has_attr(const std::string &name) const;
    bool eq_attr(const std::string &name, const gff_attr_t &val) const;
    const gff_attr_t& get_attr(const std::string &name) const;
    void set_attr_auto(const std::string &name, const std::string &val);
    void set_attr(const std::string &name, const std::string &val);
    void set_attr(const std::string &name, int64_t val);
    void set_attr(const std::string &name, double val);
    std::string str() const;
    bool empty() const;
};

struct gff_data_tmp_t {
    std::vector<gff_data_t> data;
    uint64_t linestart = 0;
    uint64_t lineend = 0;
};

class thread_pool_t
{
    struct thr_data_t {
        std::string data;
        uint64_t lnum;
    };

    std::mutex _data_mtx;
    std::mutex _tmpdata_mtx;
    int _thrmax;
    std::string *_error;
    bool _attronlystr;
    std::vector<thr_data_t> _tmpdata;
    std::vector<gff_data_tmp_t> _data;
    std::deque<std::thread> _tpool;
    static void thrproc(
        std::vector<thr_data_t> &&data, std::mutex *mtx,
        std::vector<gff_data_tmp_t> *out, std::string *error,
        bool attrval_onlystr
        );
public:
    thread_pool_t(int threads, std::string *error, bool attronlystr)
        : _thrmax(threads), _error(error), _attronlystr(attronlystr) {}
    void push(const std::string &line, uint64_t linenum);
    void flush();
    bool empty();
    std::vector<gff_data_tmp_t> data();
    void makeproc(std::vector<thr_data_t> &&data, bool sync = false);
};

class gff_parser_t
{
    int _threads;
    std::unique_ptr<thread_pool_t> _thrpool;
    int _gff_version;
    uint64_t _linenum;
    bool _onlystrval;
    std::string _error;
    std::vector<gff_data_t> _data;
    std::unordered_map<std::string, std::vector<gff_data_t*> > _data_by_seqid;
    std::unordered_map<std::string, std::vector<gff_data_t*> > _data_by_source;
    std::unordered_map<std::string, std::vector<gff_data_t*> > _data_by_type;

    enum class force_type_t{ ft_int, ft_flt, ft_str };
    std::unordered_map<std::string, force_type_t> _force_types;

public:
    static gff_data_t parse_line(const std::string &line, std::string &error, uint64_t linenum, bool onlystrinval);
    explicit gff_parser_t(int threads = 0, bool attrval_only_string = false)
        : _threads(threads)
        , _gff_version(0), _linenum(0), _onlystrval(attrval_only_string) {
        if(threads > 0) _thrpool.reset(new thread_pool_t(threads, &_error, attrval_only_string));
    }
    std::string dump(const std::string &prefix) const;
    gff_parser_t& operator<<(const std::string &line);
    bool has_error() const;
    std::string error() const;
    bool empty() const;
    uint64_t linenum() const;
    void flush();
    std::size_t size() const;

    void setattr_force_str(const std::string &field);
    void setattr_force_int(const std::string &field);
    void setattr_force_flt(const std::string &field);

    std::vector<gff_data_t*> get_by_pos(const gff_postition_t &position);
    std::vector<gff_data_t*> get_by_type(const std::string &type);
    std::vector<gff_data_t*> get_by_attr(const std::vector<gff_attribute_t> &attr);
    std::vector<gff_data_t*> get_by(
        const std::string &type,
        const std::vector<gff_attribute_t> &attr,
        const gff_postition_t &position
        );
    std::vector<gff_data_t*> get_by(
        const std::string &type,
        const std::vector<gff_attribute_t> &attr
        );
    std::vector<gff_data_t*> get_by(
        const std::string &type,
        const gff_postition_t &position
        );
    std::vector<gff_data_t*> get_by(
        const std::vector<gff_attribute_t> &attr,
        const gff_postition_t &position
        );

};

namespace utils {
bool check_no_data(const std::string &str);
std::string join_strmap(const std::unordered_map<std::string, std::string> &strmap);
std::vector<std::string> get_fields(const std::string &line, char delim, bool csv_string_format);
std::string trim(const std::string &str);
}
}

#endif // GFFPARSER_H
