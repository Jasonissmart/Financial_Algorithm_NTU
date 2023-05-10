#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
using namespace std;

class OptionCalculator {
protected:
    double S;
    double K;
    double time;
    double r;
    double sigma;

    double calculate_d1() {
        double d1 = (log(S / K) + (r + 0.5 * pow(sigma, 2)) * time) / (sigma * sqrt(time));
        return d1;
    }

    double calculate_d1(
        const double& S,
        const double& K,
        const double& time,
        const double& r,
        const double& sigma) {
        double d1 = (log(S / K) + (r + 0.5 * pow(sigma, 2)) * time) / (sigma * sqrt(time));
        return d1;
    }

    double calculate_d2() {
        double d1 = calculate_d1();
        double d2 = d1 - sigma * sqrt(time);
        return d2;
    }

    double calculate_d2(
        const double& S,
        const double& K,
        const double& time,
        const double& r,
        const double& sigma) {
        double d1 = calculate_d1(S, K, time, r, sigma);
        double d2 = d1 - sigma * sqrt(time);
        return d2;
    }

    double norm_CDF(double x) {
        return (1 + erf(x / sqrt(2))) / 2;
    }

    double norm_df(double x) {
        return exp(-0.5 * pow(x, 2)) / sqrt(2 * M_PI);
    }

public:
    OptionCalculator(const double& stock_price, const double& strike_price, const double& time_to_expiry, const double& risk_free_rate, const double& volatility) {
        this->S = stock_price;
        this->K = strike_price;
        this->time = time_to_expiry;
        this->r = risk_free_rate;
        this->sigma = volatility;
    }

    double calculate_gamma() {
        double d1 = calculate_d1();
        double gamma = norm_df(d1) / (S * sigma * sqrt(time));
        return gamma;
    }

    double calculate_vega() {
        double d1 = calculate_d1();
        double vega = S * norm_df(d1) * sqrt(time);
        return vega;
    }
};

class CallOptionCalculator : public OptionCalculator {
public:
    CallOptionCalculator(const double& stock_price, const double& strike_price, const double& time_to_expiry, const double& risk_free_rate, const double& volatility) : OptionCalculator(stock_price, strike_price, time_to_expiry, risk_free_rate, volatility)
    {
    }

    double calculate_call_price() {
        double d1 = calculate_d1();
        double d2 = calculate_d2();
        double call_price = (S * norm_CDF(d1) - K * exp(-r * time) * norm_CDF(d2));
        return call_price;
    }

    double calculate_call_price(const double& S, const double& K, const double& time, const double& r, const double& sigma) {
        double d1 = calculate_d1(S, K, time, r, sigma);
        double d2 = calculate_d2(S, K, time, r, sigma);
        double call_price = (S * norm_CDF(d1) - K * exp(-r * time) * norm_CDF(d2));
        return call_price;
    }

    double calculate_delta() {
        double d1 = calculate_d1();
        double delta = norm_CDF(d1);
        return delta;
    }

    double calculate_theta() {
        double d1 = calculate_d1();
        double d2 = calculate_d2();
        double theta = (-S * norm_df(d1) * sigma / (2 * sqrt(time)) - r * K * exp(-r * time) * norm_CDF(d2));
        return theta;
    }

    double calculate_rho() {
        double d2 = calculate_d2();
        double rho = K * time * exp(-r * time) * norm_CDF(d2);
        return rho;
    }

    void calculate_all_greeks() {
        cout << "Delta: " << this->calculate_delta() << endl;
        cout << "Gamma: " << this->calculate_gamma() << endl;
        cout << "Vega: " << this->calculate_vega() << endl;
        cout << "Theta: " << this->calculate_theta() << endl;
        cout << "Rho: " << this->calculate_rho() << endl;
    }

    // Reference: Financial Recipes 2014
    double calculate_implied_volatility(const double& option_price) {
        double vol = this->calculate_implied_volatility_newton(option_price);
        cout << "Implied Volatility:" << vol << endl;
    }

    double calculate_implied_volatility_newton(const double& option_price) {
        // check for arbitrage violations. Option price is too low if this happens
        if (option_price < 0.99 * (S - K * exp(-time * r))) {
            return 0.0;
        }
        const int MAX_ITERATION = 100;
        const double ACCURACY = 1.0e-5;
        double t_sqrt = sqrt(time);

        double newVol = (option_price / S) / (0.398 * t_sqrt); // find initial value
        for(int i = 0; i < MAX_ITERATION; i++) {
            double price = this->calculate_call_price(S, K, time, r, newVol);
            double diff = option_price - price;
            if (fabs(diff) < ACCURACY) 
                return newVol;
            double d1 = calculate_d1(S, K, time, r, newVol);
            double vega = S * t_sqrt * norm_df(d1);
            newVol = newVol + diff / vega;
        }
        return -99e10; // Something screwy happened, should throw exception
    }
};
class PutOptionCalculator : public OptionCalculator {
public:
    PutOptionCalculator(const double& stock_price, const double& strike_price, const double& time_to_expiry, const double& risk_free_rate, const double& volatility) : OptionCalculator(stock_price, strike_price, time_to_expiry, risk_free_rate, volatility)
    {
    }

    double calculate_put_price() {
        double d1 = calculate_d1();
        double d2 = calculate_d2();
        double put_price = (K * exp(-r * time) * norm_CDF(-d2) - S * norm_CDF(-d1));
        return put_price;
    }

    double calculate_delta() {
        double d1 = calculate_d1();
        double delta = norm_CDF(d1) - 1;
        return delta;
    }

    double calculate_theta() {
        double d1 = calculate_d1();
        double d2 = calculate_d2();
        double theta = (-S * norm_df(-d1) * sigma / (2 * sqrt(time)) + r * K * exp(-r * time) * norm_CDF(-d2));
        return theta;
    }

    double calculate_rho() {
        double d2 = calculate_d2();
        double rho = -K * time * exp(-r * time) * norm_CDF(-d2);
        return rho;
    }

    void calculate_all_greeks() {
        cout << "Delta: " << this->calculate_delta() << endl;
        cout << "Gamma: " << this->calculate_gamma() << endl;
        cout << "Vega: " << this->calculate_vega() << endl;
        cout << "Theta: " << this->calculate_theta() << endl;
        cout << "Rho: " << this->calculate_rho() << endl;
    }
};

int main() {
    cout << "==========================================================================" << endl;
    // Calculate black-Scholes call option price, greeks and implied volatility given market premium.
    // Inputs: stock_price, strike_price, time_to_expiry, risk_free_rate, volatility
    CallOptionCalculator callOptionCalculator(50, 50, 0.5, 0.1,  0.1);
    cout << "Call Price: " << callOptionCalculator.calculate_call_price() << endl;
    callOptionCalculator.calculate_all_greeks();
    callOptionCalculator.calculate_implied_volatility(2.5);  // Inputs: Premium

    cout << "==========================================================================" << endl;
    // Calculate black-Scholes call option price and greeks.
    // Inputs: stock_price, strike_price, time_to_expiry, risk_free_rate, volatility
    PutOptionCalculator putOptionCalculator(100, 100, 1, 0.1, 0.1);
    cout << "Put Price: " << putOptionCalculator.calculate_put_price() << endl;
    putOptionCalculator.calculate_all_greeks();
    cout << "==========================================================================" << endl;

    return 0;
}