#include <bits/stdc++.h>

//  g++ -O3 main.cpp -o main.exe

using namespace std;

using u64 = uint64_t;

struct Factor {
    u64 prime;
    u64 power;
};

struct Number {
    vector<Factor> factorization;
    Number() {
        factorization.resize(0);
    }
    Number(u64 n) {
        //factors the given number
        factorization.resize(0);
        u64 remain = n;
        for (u64 i{2}; i <= static_cast<u64>(sqrt(remain)); ++i) {
            if (remain == 1) break;
            if (remain % i == 0) {
                factorization.push_back(Factor{i, 1});
                remain /= i;
                while (remain % i == 0) {
                    ++factorization.back().power;
                    remain /= i;
                }
            }
        }
        if (remain > 1) {
            if (!factorization.empty() && factorization.back().prime == remain) {
                ++factorization.back().power;
            } else {
                factorization.push_back(Factor{remain, 1});
            }
        }
    }
    Number factor_begin() {
        Number factor(*this);
        //sets factor to 1 initially
        for (u64 i{}; i < factor.factorization.size(); ++i) {
            Factor& f = factor.factorization[i];
            f.power = 0;
        }
        return factor;
    }
};

ostream& operator<<(ostream& out, Number n) {
    //prints prime factorization
    for (u64 i{}; i < n.factorization.size(); ++i) {
        Factor& f = n.factorization[i];
        out << f.prime << '^' << f.power;
        if (i != n.factorization.size() - 1) out << " * ";
    }
    return out;
}

u64 value(Number n) {
    //the number being represented
    u64 res = 1;
    for (u64 i{}; i < n.factorization.size(); ++i) {
        Factor& f = n.factorization[i];
        if (f.power == 0) continue;
        res *= pow(f.prime, f.power);
    }
    return res;
}

u64 phi(Number n) {
    //computes euler totient function using factorization
    u64 res = 1;
    for (u64 i{}; i < n.factorization.size(); ++i) {
        Factor& f = n.factorization[i];
        if (f.power == 0) continue;
        res *= (f.prime - 1) * pow(f.prime, f.power - 1);
    }
    return res;
}

u64 sigma(Number n) {
    //computes sum of divisors function using factorization
    u64 res = 1;
    for (u64 i{}; i < n.factorization.size(); ++i) {
        Factor& f = n.factorization[i];
        if (f.power == 0) continue;
        res *= (pow(f.prime, f.power + 1) - 1) / (f.prime - 1);
    }
    return res;
}

bool next_factor(Number& factor, Number n) {
    assert(factor.factorization.size() == n.factorization.size());
    int i = 0;
    while (true) {
        factor.factorization[i].power++;
        if (factor.factorization[i].power > n.factorization[i].power) {
            factor.factorization[i].power = 0;
            ++i;
            if (i == factor.factorization.size()) {
                return true;
            }
        } else {
            break;
        }
    }
    return false;
}

u64 divisor_sum(u64 input) {
    Number n(input);
    Number factor(n);
    //sets factor to 1 initially
    for (u64 i{}; i < factor.factorization.size(); ++i) {
        Factor& f = factor.factorization[i];
        f.power = 0;
    }
    //sums the expression
    u64 sum{0};
    while (true) {
        //cout << "factor: " << factor << '\n';
        u64 term = value(factor);
        sum += term;
        bool done = next_factor(factor, n);
        if (done) break;
    }
    return sum;
}

bool find_hybrid_2np(u64 input, std::ofstream& fout) {
    Number n(input);
    u64 sigma_n = sigma(n);
    for(Number d1 = n.factor_begin();;) {
        //cout << "factor: " << factor << '\n';
        //process stuff
        for (Number d2 = n.factor_begin();;) {
            if (value(d1) != value(d2)) {
                u64 vd1 = value(d1), vd2 = value(d2);
                if (sigma_n == 2 * input + vd1 + vd2) {
                    fout << "sigma(" << input << ") = " << sigma_n << " = 2 * " << input << " + " << vd1 << " + " << vd2 << std::endl;
                    return true;
                }
                if (sigma_n == 2 * input + vd1 - vd2) {
                    fout << "sigma(" << input << ") = " << sigma_n << " = 2 * " << input << " + " << vd1 << " - " << vd2 << std::endl;
                    return true;
                }
                if (sigma_n == 2 * input - vd1 - vd2) {
                    
                    fout << "sigma(" << input << ") = " << sigma_n << " = 2 * " << input << " - " << vd1 << " - " << vd2 << std::endl;
                    return true;
                }
            }
            bool done = next_factor(d2, n);
            if (done) break;
        }
        bool done = next_factor(d1, n);
        if (done) break;
    }
    //cout << input << " not hybrid 2-near perfect\n";
    return false;
}

int main() {
    std::ofstream fout("output.txt");
    //u64 input;
    //cin >> input;
    fout << setprecision(12);
    for (int i{2}; i < 10000000; ++i) {
        Number n(i);
        //if ((n.factorization.size() >= 2 || n.factorization[0].power >= 2)) {
        if (n.factorization.size() == 2 && n.factorization[0].prime == 2) {
            bool twonp = find_hybrid_2np(i, fout);
            if (twonp) {
                fout << "----------------------------------------------------->" << n << std::endl;
            }
        }
    }
}
