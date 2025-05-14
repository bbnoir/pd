#include "ObjectiveFunction.h"

#include "cstdio"
#include "cmath"

using namespace std;

ExampleFunction::ExampleFunction(Placement &placement) : BaseFunction(1), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &ExampleFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &ExampleFunction::Backward()
{
    // Compute the gradient of the function
    grad_[0].x = 6. * input_[0].x + 2. * input_[0].y;
    grad_[0].y = 2. * input_[0].x + 4. * input_[0].y;
    return grad_;
}

WAWirelength::WAWirelength(Placement &placement) : BaseFunction(placement.numModules()), placement_(placement)
{
    gamma_ = (placement_.boundryRight() - placement_.boundryLeft()) * 0.01;
    const size_t num_nets = placement_.numNets();
    x_max_.resize(num_nets);
    y_max_.resize(num_nets);
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
    input_ = input;
    value_ = 0.;
    for (size_t i = 0, end_i = placement_.numNets(); i < end_i; ++i) {
        Net &net = placement_.net(i);

        // maximum and minimum x and y
        double x_max = placement_.boundryLeft();
        double y_max = placement_.boundryBottom();
        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Pin &pin = net.pin(j);
            Point2<double> &module_pos = input_[pin.moduleId()];
            x_max = max(x_max, module_pos.x + pin.xOffset());
            y_max = max(y_max, module_pos.y + pin.yOffset());
        }
        // cache
        x_max_[i] = x_max;
        y_max_[i] = y_max;

        // sum of exponential
        double x_max_sum_w_exp = 0., x_max_sum_exp = 0.;
        double x_min_sum_w_exp = 0., x_min_sum_exp = 0.;
        double y_max_sum_w_exp = 0., y_max_sum_exp = 0.;
        double y_min_sum_w_exp = 0., y_min_sum_exp = 0.;
        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Pin &pin = net.pin(j);
            Point2<double> &module_pos = input_[pin.moduleId()];
            const double x = module_pos.x + pin.xOffset();
            const double y = module_pos.y + pin.yOffset();
            const double exp_x_max = exp((x - x_max) / gamma_);
            const double exp_x_min = exp((x_max - x) / gamma_);
            const double exp_y_max = exp((y - y_max) / gamma_);
            const double exp_y_min = exp((y_max - y) / gamma_);
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
        const double y_max = y_max_[i];
        const double x_max_sum_w_exp = x_max_sum_w_exp_[i];
        const double x_min_sum_w_exp = x_min_sum_w_exp_[i];
        const double y_max_sum_w_exp = y_max_sum_w_exp_[i];
        const double y_min_sum_w_exp = y_min_sum_w_exp_[i];
        const double x_max_sum_exp = x_max_sum_exp_[i];
        const double x_min_sum_exp = x_min_sum_exp_[i];
        const double y_max_sum_exp = y_max_sum_exp_[i];
        const double y_min_sum_exp = y_min_sum_exp_[i];

        for (size_t j = 0, end_j = net.numPins(); j < end_j; ++j) {
            Pin &pin = net.pin(j);
            Point2<double> &module_pos = input_[pin.moduleId()];
            const double x = module_pos.x + pin.xOffset();
            const double y = module_pos.y + pin.yOffset();
            const double exp_x_max = exp((x - x_max) / gamma_);
            const double exp_x_min = exp((x_max - x) / gamma_);
            const double exp_y_max = exp((y - y_max) / gamma_);
            const double exp_y_min = exp((y_max - y) / gamma_);
            const double grad_x_max = ((1 + x / gamma_) * exp_x_max * x_max_sum_exp - x_max_sum_w_exp / gamma_ * exp_x_max) / x_max_sum_exp / x_max_sum_exp;
            const double grad_x_min = ((1 - x / gamma_) * exp_x_min * x_min_sum_exp + x_min_sum_w_exp / gamma_ * exp_x_min) / x_min_sum_exp / x_min_sum_exp;
            const double grad_y_max = ((1 + y / gamma_) * exp_y_max * y_max_sum_exp - y_max_sum_w_exp / gamma_ * exp_y_max) / y_max_sum_exp / y_max_sum_exp;
            const double grad_y_min = ((1 - y / gamma_) * exp_y_min * y_min_sum_exp + y_min_sum_w_exp / gamma_ * exp_y_min) / y_min_sum_exp / y_min_sum_exp;
            auto &grad = grad_[pin.moduleId()];
            grad.x += grad_x_max - grad_x_min;
            grad.y += grad_y_max - grad_y_min;
        }
    }
    return grad_;
}

Density::Density(Placement &placement) : BaseFunction(placement.numModules()), placement_(placement)
{
}

const double &Density::operator()(const std::vector<Point2<double>> &input)
{
    input_ = input;
    value_ = 0.;
    // TODO: Implement the density function
    return value_;
}

const std::vector<Point2<double>> &Density::Backward()
{
    // TODO: Implement the backward pass of the density function
    // reset gradient
    for (auto &grad : grad_) {
        grad.x = 0.;
        grad.y = 0.;
    }
    return grad_;
}

ObjectiveFunction::ObjectiveFunction(Placement &placement) : BaseFunction(placement.numModules()), wirelength_(placement), density_(placement)
{
    lambda_ = 0.1;
}

const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input)
{
    value_ = wirelength_(input) + lambda_ * density_(input);
    return value_;
}

const std::vector<Point2<double>> &ObjectiveFunction::Backward()
{
    wirelength_.Backward();
    density_.Backward();
    for (size_t i = 0, end_i = grad_.size(); i < end_i; ++i) {
        grad_[i] = wirelength_.grad()[i] + lambda_ * density_.grad()[i];
    }
    return grad_;
}


