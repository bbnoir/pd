#define _GLIBCXX_USE_CXX11_ABI 0  // Align the ABI version to avoid compatibility issues with `Placment.h`
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <vector>

#include "Placement.h"
#include "Point.h"
#include "Rectangle.h"

/**
 * @brief Base class for objective functions
 */
class BaseFunction {
   public:
    /////////////////////////////////
    // Conssutructors
    /////////////////////////////////

    BaseFunction(const size_t &input_size) : grad_(input_size) {}

    /////////////////////////////////
    // Accessors
    /////////////////////////////////

    const std::vector<Point2<double>> &grad() const { return grad_; }
    const double &value() const { return value_; }

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    // Forward pass, compute the value of the function
    virtual const double &operator()(const std::vector<Point2<double>> &input) = 0;

    // Backward pass, compute the gradient of the function
    virtual const std::vector<Point2<double>> &Backward() = 0;

    virtual double bin_width() const { return 0.; }

   protected:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> grad_;  // Gradient of the function
    double value_;                      // Value of the function
};

class WAWirelength : public BaseFunction {
public:
    WAWirelength(Placement &placement, std::vector<Point2<double>> &input);
    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;
private:
    Placement &placement_;
    std::vector<Point2<double>> &input_;
    double gamma_;

    std::vector<double> x_max_, y_max_, x_min_, y_min_;
    std::vector<double> x_max_sum_w_exp_, x_min_sum_w_exp_;
    std::vector<double> y_max_sum_w_exp_, y_min_sum_w_exp_;
    std::vector<double> x_max_sum_exp_, x_min_sum_exp_;
    std::vector<double> y_max_sum_exp_, y_min_sum_exp_;
};

/**
 * @brief Density function
 */
class Density : public BaseFunction {
    // TODO: Implement the density function, add necessary data members for caching
   public:
    /////////////////////////////////
    // Methods
    /////////////////////////////////
    Density(Placement &placement, std::vector<Point2<double>> &input);

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

    void resize_bin(int grid_num);

    double bin_width() const { return bin_width_; }
    double overflow_ratio() const { return overflow_ratio_; }

   private:
    Placement &placement_;
    std::vector<Point2<double>> &input_;
    size_t num_bins_x_;
    size_t num_bins_y_;
    double bin_width_;
    double bin_height_;
    double bin_area_;
    double t_density_;
    int bin_range_;
    double total_moveable_area_;
    double overflow_ratio_;
    std::vector<std::vector<double>> Mb_, Db_;

};

/**
 * @brief Objective function for global placement
 */
class ObjectiveFunction : public BaseFunction {
    // TODO: Implement the objective function for global placement, add necessary data
    // members for caching
    //
    // Hint: The objetive function of global placement is as follows:
    //       f(t) = wirelength(t) + lambda * density(t),
    // where t is the positions of the modules, and lambda is the penalty weight.
    // You may need an interface to update the penalty weight (lambda) dynamically.
   public:
    /////////////////////////////////
    // Methods
    /////////////////////////////////
    ObjectiveFunction(Placement &placement, std::vector<Point2<double>> &input);

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

    void resize_bin(int grid_num) { density_.resize_bin(grid_num); }
    double bin_width() const { return density_.bin_width(); }
    double overflow_ratio() const { return density_.overflow_ratio(); } 
    double lambda() const { return lambda_; }
    double wl_value() const { return wirelength_.value(); }
    double d_value() const { return density_.value(); }

    void set_init_lambda();
    inline void scale_lambda(double scale) { lambda_ *= scale; }
    inline void set_lambda(double lambda) { lambda_ = lambda; }

   private:
    WAWirelength wirelength_;
    Density density_;
    std::vector<Point2<double>> &input_;
    double lambda_;
};

#endif  // OBJECTIVEFUNCTION_H
