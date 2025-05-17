#include "ObjectiveFunction.h"

#include "cstdio"
#include "cmath"
#include "Rectangle.h"

using namespace std;

WAWirelength::WAWirelength(Placement &placement, std::vector<Point2<double>> &input)
    : BaseFunction(placement.numModules()), placement_(placement), input_(input)
{
    gamma_ = max(placement_.boundryRight() - placement_.boundryLeft(),
                  placement_.boundryTop() - placement_.boundryBottom()) * 0.1;
    const size_t num_nets = placement_.numNets();
    x_max_.resize(num_nets);
    y_max_.resize(num_nets);
    x_min_.resize(num_nets);
    y_min_.resize(num_nets);
    x_max_sum_w_exp_.resize(num_nets);
    x_min_sum_w_exp_.resize(num_nets);
    y_max_sum_w_exp_.resize(num_nets);
    y_min_sum_w_exp_.resize(num_nets);
    x_max_sum_exp_.resize(num_nets);
    x_min_sum_exp_.resize(num_nets);
    y_max_sum_exp_.resize(num_nets);
    y_min_sum_exp_.resize(num_nets);
}

const double &WAWirelength::operator()(const std::vector<Point2<double>> &input) 
{
    value_ = 0.;
    for (size_t i = 0, end_i = placement_.numNets(); i < end_i; ++i) {
        Net &net = placement_.net(i);

        // maximum and minimum x and y
        double x_max = placement_.boundryLeft();
        double x_min = placement_.boundryRight();
        double y_max = placement_.boundryBottom();
        double y_min = placement_.boundryTop();
        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Point2<double> module_pos = input_[net.pin(j).moduleId()];
            x_max = max(x_max, module_pos.x);
            x_min = min(x_min, module_pos.x);
            y_max = max(y_max, module_pos.y);
            y_min = min(y_min, module_pos.y);
        }
        // cache
        x_max_[i] = x_max;
        x_min_[i] = x_min;
        y_max_[i] = y_max;
        y_min_[i] = y_min;

        // sum of exponential
        double x_max_sum_w_exp = 0., x_max_sum_exp = 0.;
        double x_min_sum_w_exp = 0., x_min_sum_exp = 0.;
        double y_max_sum_w_exp = 0., y_max_sum_exp = 0.;
        double y_min_sum_w_exp = 0., y_min_sum_exp = 0.;
        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Point2<double> &module_pos = input_[net.pin(j).moduleId()];
            const double x = module_pos.x;
            const double y = module_pos.y;
            const double exp_x_max = exp((x - x_max) / gamma_);
            const double exp_x_min = exp((x_min - x) / gamma_);
            const double exp_y_max = exp((y - y_max) / gamma_);
            const double exp_y_min = exp((y_min - y) / gamma_);
            x_max_sum_w_exp += exp_x_max * x;
            x_min_sum_w_exp += exp_x_min * x;
            y_max_sum_w_exp += exp_y_max * y;
            y_min_sum_w_exp += exp_y_min * y;
            x_max_sum_exp += exp_x_max;
            x_min_sum_exp += exp_x_min;
            y_max_sum_exp += exp_y_max;
            y_min_sum_exp += exp_y_min;
        }
        // cache
        x_max_sum_w_exp_[i] = x_max_sum_w_exp;
        x_min_sum_w_exp_[i] = x_min_sum_w_exp;
        y_max_sum_w_exp_[i] = y_max_sum_w_exp;
        y_min_sum_w_exp_[i] = y_min_sum_w_exp;
        x_max_sum_exp_[i] = x_max_sum_exp;
        x_min_sum_exp_[i] = x_min_sum_exp;
        y_max_sum_exp_[i] = y_max_sum_exp;
        y_min_sum_exp_[i] = y_min_sum_exp;

        // wirelength
        value_ += x_max_sum_w_exp / x_max_sum_exp - x_min_sum_w_exp / x_min_sum_exp +
                  y_max_sum_w_exp / y_max_sum_exp - y_min_sum_w_exp / y_min_sum_exp;
    }

    return value_;
}

const std::vector<Point2<double>> &WAWirelength::Backward()
{
    // reset gradient
    for (auto &grad : grad_) {
        grad.x = 0.;
        grad.y = 0.;
    }
    for (size_t i = 0, end_i = placement_.numNets(); i < end_i; ++i) {
        Net &net = placement_.net(i);
        const double x_max = x_max_[i];
        const double x_min = x_min_[i];
        const double y_max = y_max_[i];
        const double y_min = y_min_[i];
        const double x_max_sum_w_exp = x_max_sum_w_exp_[i];
        const double x_min_sum_w_exp = x_min_sum_w_exp_[i];
        const double y_max_sum_w_exp = y_max_sum_w_exp_[i];
        const double y_min_sum_w_exp = y_min_sum_w_exp_[i];
        const double x_max_sum_exp = x_max_sum_exp_[i];
        const double x_min_sum_exp = x_min_sum_exp_[i];
        const double y_max_sum_exp = y_max_sum_exp_[i];
        const double y_min_sum_exp = y_min_sum_exp_[i];

        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Point2<double> &module_pos = input_[net.pin(j).moduleId()];
            const double x = module_pos.x;
            const double y = module_pos.y;
            const double exp_x_max = exp((x - x_max) / gamma_);
            const double exp_x_min = exp((x_min - x) / gamma_);
            const double exp_y_max = exp((y - y_max) / gamma_);
            const double exp_y_min = exp((y_min - y) / gamma_);
            const double grad_x_max = ((1 + x / gamma_) * exp_x_max * x_max_sum_exp - x_max_sum_w_exp / gamma_ * exp_x_max) / x_max_sum_exp / x_max_sum_exp;
            const double grad_x_min = ((1 - x / gamma_) * exp_x_min * x_min_sum_exp + x_min_sum_w_exp / gamma_ * exp_x_min) / x_min_sum_exp / x_min_sum_exp;
            const double grad_y_max = ((1 + y / gamma_) * exp_y_max * y_max_sum_exp - y_max_sum_w_exp / gamma_ * exp_y_max) / y_max_sum_exp / y_max_sum_exp;
            const double grad_y_min = ((1 - y / gamma_) * exp_y_min * y_min_sum_exp + y_min_sum_w_exp / gamma_ * exp_y_min) / y_min_sum_exp / y_min_sum_exp;
            auto &grad = grad_[net.pin(j).moduleId()];
            grad.x += grad_x_max - grad_x_min;
            grad.y += grad_y_max - grad_y_min;
        }
    }
    return grad_;
}

Density::Density(Placement &placement, std::vector<Point2<double>> &input)
    : BaseFunction(placement.numModules()), placement_(placement), input_(input)
{
    // initialize bins
    bin_range_ = 2;
    t_density_ = 0.9;
    const int grid_num = 8;
    this->resize_bin(grid_num);
    total_moveable_area_ = 0.;
    for (size_t i = 0, end_i = placement_.numModules(); i < end_i; ++i) {
        Module &module = placement_.module(i);
        if (!module.isFixed()) {
            total_moveable_area_ += module.width() * module.height();
        }
    }
}

void Density::resize_bin(int grid_num)
{
    overflow_ratio_ = 0.;
    num_bins_x_ = num_bins_y_ = grid_num;
    bin_width_ = (placement_.boundryRight() - placement_.boundryLeft()) / grid_num;
    bin_height_ = (placement_.boundryTop() - placement_.boundryBottom()) / grid_num;
    bin_area_ = bin_width_ * bin_height_;
    Mb_.clear();
    Db_.clear();
    Mb_.resize(num_bins_y_, vector<double>(num_bins_x_, bin_area_));
    Db_.resize(num_bins_y_, vector<double>(num_bins_x_, 0.));
    for (size_t i = 0, end_i = placement_.numModules(); i < end_i; ++i) {
        Module &module = placement_.module(i);
        if (module.isFixed()) {
            Rectangle module_rect = module.rectangle();
            const double x_min = module_rect.left();
            const double y_min = module_rect.bottom();
            const double x_max = module_rect.right();
            const double y_max = module_rect.top();
            int bin_x_min = (x_min - placement_.boundryLeft()) / bin_width_;
            bin_x_min = max(bin_x_min, 0);
            int bin_y_min = (y_min - placement_.boundryBottom()) / bin_height_;
            bin_y_min = max(bin_y_min, 0);
            int bin_x_max = (x_max - placement_.boundryLeft()) / bin_width_;
            bin_x_max = min(bin_x_max, int(num_bins_x_ - 1));
            int bin_y_max = (y_max - placement_.boundryBottom()) / bin_height_;
            bin_y_max = min(bin_y_max, int(num_bins_y_ - 1));
            for (int j = bin_x_min; j <= bin_x_max; ++j) {
                for (int k = bin_y_min; k <= bin_y_max; ++k) {
                    const Rectangle bin_rect = Rectangle(placement_.boundryLeft() + j * bin_width_,
                                                placement_.boundryBottom() + k * bin_height_,
                                                placement_.boundryLeft() + (j + 1) * bin_width_,
                                                placement_.boundryBottom() + (k + 1) * bin_height_);
                    Mb_[j][k] -= Rectangle::overlapArea(bin_rect, module_rect);
                }
            }
        }
    }
    for (size_t i = 0, end_i = num_bins_y_; i < end_i; ++i) {
        for (size_t j = 0, end_j = num_bins_x_; j < end_j; ++j) {
            Mb_[i][j] *= t_density_;
        }
    }
}

static inline double bell_shaped(double dx, double wv, double wb, double a, double b)
{
    double px = 0;
    if (dx <= wv * 0.5 + wb) {
        px = 1. - a * dx * dx;
    } else if ((wv * 0.5 + wb) < dx && dx <= (wv * 0.5 + 2. * wb)) {
        px = b * pow(dx - wv * 0.5 - wb * 2., 2);
    }
    return px;
}

static inline double d_bell_shaped(double dx, double wv, double wb, double a, double b)
{
    const double abs_dx = abs(dx);
    double dpx = 0;
    if (abs_dx <= wv * 0.5 + wb) {
        dpx = -2. * a * dx;
    } else if ((wv * 0.5 + wb) < abs_dx && abs_dx <= (wv * 0.5 + 2. * wb)) {
        dpx = (dx > 0) ? 2. * b * (dx - wv * 0.5 - wb * 2.) : 2. * b * (dx + wv * 0.5 + wb * 2.);
    }
    return dpx;
}

const double &Density::operator()(const std::vector<Point2<double>> &input)
{
    value_ = 0.;
    // reset density
    for (size_t i = 0, end_i = num_bins_y_; i < end_i; ++i) {
        for (size_t j = 0, end_j = num_bins_x_; j < end_j; ++j) {
            Db_[i][j] = 0.;
        }
    }
    vector<vector<double>> tempDb(num_bins_y_, vector<double>(num_bins_x_, 0.));
    // calculate density
    for (size_t i = 0, end_i = placement_.numModules(); i < end_i; ++i) {
        Module &module = placement_.module(i);
        const double w_v = module.width();
        const double h_v= module.height();
        const double w_b = bin_width_;
        const double h_b = bin_height_;
        const double lx = input_[i].x;
        const double rx = input_[i].x + w_v;
        const double by = input_[i].y;
        const double ty = input_[i].y + h_v;
        int bin_lc = (lx - placement_.boundryLeft()) / bin_width_;
        bin_lc = max(bin_lc - bin_range_, 0);
        int bin_rc = (rx - placement_.boundryLeft()) / bin_width_;
        bin_rc = min(bin_rc + bin_range_, int(num_bins_x_ - 1));
        int bin_bc = (by - placement_.boundryBottom()) / bin_height_;
        bin_bc = max(bin_bc - bin_range_, 0);
        int bin_tc = (ty - placement_.boundryBottom()) / bin_height_;
        bin_tc = min(bin_tc + bin_range_, int(num_bins_y_ - 1));
        // bell-shaped
        const double a_x = 4. / (w_v + 2. * w_b) / (w_v + 4. * w_b);
        const double a_y = 4. / (h_v + 2. * h_b) / (h_v + 4. * h_b);
        const double b_x = 2. / w_b / (w_v + 4. * w_b);
        const double b_y = 2. / h_b / (h_v + 4. * h_b);
        double total_potential = 0.;
        for (int j = bin_lc; j <= bin_rc; ++j) {
            for (int k = bin_bc; k <= bin_tc; ++k) {
                const double bin_center_x = placement_.boundryLeft() + j * bin_width_ + bin_width_ / 2.;
                const double dx = abs(bin_center_x - (lx + w_v / 2.));
                const double px = bell_shaped(dx, w_v, w_b, a_x, b_x);
                const double bin_center_y = placement_.boundryBottom() + k * bin_height_ + bin_height_ / 2.;
                const double dy = abs(bin_center_y - (by + h_v / 2.));
                const double py = bell_shaped(dy, h_v, h_b, a_y, b_y);
                total_potential += px * py;
                tempDb[j][k] = px * py;
            }
        }
        // normalize
        const double module_area = w_v * h_v;
        for (int j = bin_lc; j <= bin_rc; ++j) {
            for (int k = bin_bc; k <= bin_tc; ++k) {
                Db_[j][k] += tempDb[j][k] / total_potential * module_area;
            }
        }
    }
    overflow_ratio_ = 0.;
    for (size_t i = 0, end_i = num_bins_y_; i < end_i; ++i) {
        for (size_t j = 0, end_j = num_bins_x_; j < end_j; ++j) {
            value_ += pow(Db_[i][j] - Mb_[i][j], 2);
            overflow_ratio_ += max(Db_[i][j] - Mb_[i][j], 0.);
        }
    }
    overflow_ratio_ /= total_moveable_area_;
    return value_;
}

const std::vector<Point2<double>> &Density::Backward()
{
    // reset gradient
    for (auto &grad : grad_) {
        grad.x = 0.;
        grad.y = 0.;
    }
    // calculate gradient
    for (size_t i = 0, end_i = placement_.numModules(); i < end_i; ++i) {
        Module &module = placement_.module(i);
        const double w_v = module.width();
        const double h_v= module.height();
        const double w_b = bin_width_;
        const double h_b = bin_height_;
        const double lx = input_[i].x;
        const double rx = input_[i].x + w_v;
        const double by = input_[i].y;
        const double ty = input_[i].y + h_v;
        int bin_lc = (lx - placement_.boundryLeft()) / bin_width_;
        bin_lc = max(bin_lc - bin_range_, 0);
        int bin_rc = (rx - placement_.boundryLeft()) / bin_width_;
        bin_rc = min(bin_rc + bin_range_, int(num_bins_x_ - 1));
        int bin_bc = (by - placement_.boundryBottom()) / bin_height_;
        bin_bc = max(bin_bc - bin_range_, 0);
        int bin_tc = (ty - placement_.boundryBottom()) / bin_height_;
        bin_tc = min(bin_tc + bin_range_, int(num_bins_y_ - 1));
        // bell-shaped
        const double a_x = 4. / (w_v + 2. * w_b) / (w_v + 4. * w_b);
        const double a_y = 4. / (h_v + 2. * h_b) / (h_v + 4. * h_b);
        const double b_x = 2. / w_b / (w_v + 4. * w_b);
        const double b_y = 2. / h_b / (h_v + 4. * h_b);
        for (int j = bin_lc; j <= bin_rc; ++j) {
            for (int k = bin_bc; k <= bin_tc; ++k) {
                const double bin_center_x = placement_.boundryLeft() + j * bin_width_ + bin_width_ / 2.;
                const double dx = (lx + w_v / 2.) - bin_center_x;
                const double px = bell_shaped(abs(dx), w_v, w_b, a_x, b_x);
                const double dpx_dx = d_bell_shaped(dx, w_v, w_b, a_x, b_x);
                const double bin_center_y = placement_.boundryBottom() + k * bin_height_ + bin_height_ / 2.;
                const double dy = (by + h_v / 2.) - bin_center_y;
                const double py = bell_shaped(abs(dy), h_v, h_b, a_y, b_y);
                const double dpy_dy = d_bell_shaped(dy, h_v, h_b, a_y, b_y);
                auto &grad = grad_[i];
                const double db_mb = 2. * (Db_[j][k] - Mb_[j][k]);
                grad.x += db_mb * py * dpx_dx;
                grad.y += db_mb * px * dpy_dy;
            }
        }
    }
    return grad_;
}

ObjectiveFunction::ObjectiveFunction(Placement &placement, std::vector<Point2<double>> &input)
    : BaseFunction(placement.numModules()), wirelength_(placement, input), density_(placement, input), input_(input)
{
    lambda_ = 0;
}

const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input)
{
    value_ = wirelength_(input);
    if (lambda_ != 0) { value_ += lambda_ * density_(input); }
    return value_;
}

const std::vector<Point2<double>> &ObjectiveFunction::Backward()
{
    wirelength_.Backward();
    if (lambda_ != 0) { density_.Backward(); }
    for (size_t i = 0, end_i = grad_.size(); i < end_i; ++i) {
        grad_[i] = wirelength_.grad()[i];
    }
    if (lambda_ != 0) {
        for (size_t i = 0, end_i = grad_.size(); i < end_i; ++i) {
            grad_[i] += lambda_ * density_.grad()[i];
        }
    }
    return grad_;
}

void ObjectiveFunction::set_init_lambda()
{
    value_ = wirelength_(input_) + lambda_ * density_(input_);
    wirelength_.Backward();
    density_.Backward();
    double sum_wl_grad = 0., sum_d_grad = 0.;
    for (auto &grad : wirelength_.grad()) { sum_wl_grad += Norm2(grad); }
    for (auto &grad : density_.grad()) { sum_d_grad += Norm2(grad); }
    lambda_ = sum_wl_grad / sum_d_grad;
}