#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <random>

using namespace std;

// Граница для сравнения с единицей
double eps_sr = 1e-14;

// ПУНКТ 1
// Сформируем функции генерации регулярных и адаптивных сеточных разбиений произвольного отрезка [a,b] в зависимости от числа сегментов разбиения и величины коэффициента разрядки r

// Функция для формирования равномерного разбиения

vector<double> genRegGrd(double const &a, double const &b, int n){

    if(n < 1) {
        throw invalid_argument("Количество сегментов должно быть >=1");
    }

    // Форимруем шаг сетки, где указанный n это число сегментов разбиения
    double h = (b - a) / n;

    // Выделяем память под массив узлов
    vector<double> x(n + 1);

    for (int i = 0; i <= n; i++) {

        if (i != n) {
            x[i] = a + i * h;
        }else {
            x[i] = b;
        }
    }

    return x;
}

// Опишем функцию для формирования адаптивного сеточного разбиения

vector<double> genAdpGrd(double const &a, double const &b, double r, int n){

    if(n < 1) {
        throw invalid_argument("Количество сегментов должно быть >=1");
    }

    if(fabs(r - 1.0) < eps_sr) {
        throw invalid_argument("Используйте функцию для построения равномерного разбиения");
    }

    // Для формирования адаптивного шага сетки можем воспользоваться формулой геометрической прогрессии
    // Sn=b1*(1-q^(n))/(1-q), где q=r - > b-a = b1*(1-r^(n))/(1-r) -> b1 = (b-a)*(1-r)/(1-r^(n))

    // Форимруем начальный шаг сетки, где указанный n это число сегментов разбиения
    double h = (b - a) * (1 - r) / (1 - pow(r, n));

    // Выделяем память под массив узлов
    vector<double> x(n + 1);

    x[0] = a;

    for (int i = 1; i <= n; i++) {

        if (i != n) {
            x[i] = x[i - 1] + h * pow(r, i - 1);
        }else {
            x[i] = b;
        }
    }
    return x;
}

// ПУНКТ 2
// Опишем класс реализующий интерфейс кубического интерполяционного сплайна

class Cubic_Int_Spl {

    private:
        vector<double> x;
        vector<double> a, b, c, d;

    public:

        void update_spl(const vector<double> &x_upd, const vector<double> &f_d){

            x.clear();

            for(auto &elem :x_upd) x.push_back(elem);

            // определяем число отрезков разбиения
            int n = static_cast<int>(x.size()) - 1;

            double h_c, h_n;

            // изменяем размеры векторов коэффициентов

            a.resize(n);
            b.resize(n);
            c.resize(n);
            d.resize(n);

            // сформируем вектор правой части СЛАУ

            vector<double>f(n - 1);

            // вычислим соответсвующие коэффициенты по формулам из учебного пособия
            for (int i = 0; i < n - 1; i++) {

                h_c = x[i + 1] - x[i];
                h_n = x[i + 2] - x[i + 1];
                // формируем главную диагональ
                b[i] = 2 * (h_c + h_n);
                // формируем нижнюю диагональ
                a[i + 1] = h_c;
                // формируем верхнюю диагональ
                d[i] = h_n;
                // правая часть СЛАУ
                f[i] = 3.0 * ((f_d[i + 2] - f_d[i + 1]) / h_n - (f_d[i + 1] - f_d[i]) / h_c);

            }

            // реализуем прямой и обратный ход методом прогонки

            for(int j = 1; j < n - 1; j++) {
                b[j] -= a[j] / b[j - 1] * d[j - 1];
                f[j] -= a[j] / b[j - 1] * f[j - 1];
            }

            c[n - 1] = f[n - 2] / b[n - 2];

            for (int j = n - 2; j > 0; j--) {
                c[j] = (f[j - 1] - c[j + 1] * d[j - 1]) / b[j - 1];
            }

            // формируем условия нулевой кривизны

            c[0] = 0.0;

            // формируем коэффиценты сплайна

            for(int i = 0; i < n - 1; i++) {
                h_c = x[i + 1] - x[i];
                a[i] = f_d[i];
                b[i] = (f_d[i + 1] - f_d[i]) / h_c - (c[i + 1] + 2.0 * c[i]) * h_c / 3.0;
                d[i] = (c[i + 1] - c[i]) / h_c / 3.0;
            }


            // последний сегмент обрабатываем отдельно
            h_c = x[n] - x[n - 1];
            a[n - 1] = f_d[n - 1];
            b[n - 1] = (f_d[n] - f_d[n - 1]) / h_c - 2.0 * c[n - 1] * h_c / 3.0;
            d[n - 1] = -c[n - 1] / h_c / 3.0;

        }

        vector<double> get_spl_and_two_dev(const double &x_in){

            // определяем число отрезков разбиения
            int n = static_cast<int>(x.size()) - 1;

            vector<double> res(3);

            for(int i = 0; i < n; i++) {
                if((x_in > x[i] && x_in < x[i + 1]) || fabs(x_in - x[i]) < eps_sr || fabs(x_in - x[i + 1]) < eps_sr) {
                    double local = x_in - x[i];
                    res[0] = a[i] + b[i] * local + c[i] * pow(local, 2) + d[i] * pow(local, 3);
                    res[1] = b[i] + 2.0 * c[i] * local + 3.0 * d[i] * pow(local, 2);
                    res[2] = 2.0 * c[i] + 6.0 * d[i] * local;

                    return res;
                }

            }
            throw invalid_argument("Указано не корректное число");
        }

        // функция вычисления погрешности вычислений
        vector<double> err_of_aprx(){

            // в рамках данного ПЗ будем вычислять погрешность используя норму пространства C[a,b]

            int N = 100000;

            double err;
            double err_d;
            double err_d2;
            double max_err_s = 0;
            double max_err_s_dev = 0;
            double max_err_s_dev2 = 0;

            vector<double>grd = genRegGrd(x[0], x[x.size() - 1], N);

            for(int i = 0; i < N; i++) {

                vector<double>res = get_spl_and_two_dev(grd[i]);

                err = fabs(res[0] - sin(grd[i]));
                err_d = fabs(res[1] - cos(grd[i]));
                err_d2 = fabs(res[2] - (-sin(grd[i])));

                if(max_err_s < err) {
                    max_err_s = err;
                }
                if(max_err_s_dev < err_d) {
                    max_err_s_dev = err_d;
                }
                if(max_err_s_dev2 < err_d2) {
                    max_err_s_dev2 = err_d2;
                }

            }

            vector<double>res_err = {
                max_err_s, max_err_s_dev, max_err_s_dev2
            };

            return res_err;
        }

        vector<double>genRdGrd_reg(){

            random_device rd;

            mt19937 gen(rd());

            uniform_real_distribution<double> dist(0.1, 0.9);

            vector<double>x_rd(x.size() - 1);

            double h = x[1] - x[0];

            for(int i = 0; i < x_rd.size(); i++) {
                double local = h * dist(gen);
                x_rd[i] = x[0] + h * i + local;
            }

            return x_rd;
        }

};



int main(int argc, const char *argv[]) {

    // ПУНКТ 3

    double a, b;
    int n;

    cout << "Введите границы отрезка a и b: ";
    cin >> a >> b;
    cout << "Введите число сегментов разбиения: ";
    cin >> n;

    // поскольку далее нам будет необходимо построить не менее 10 точек не совпадающих с узловыми введем ограничение

    if (n <= 10) {
        throw invalid_argument("Количество сегментов должно быть > 10");
    }

    // в качестве непрерывной неполиномиальной функции будем использовать sin(x);

    // формируем равномерное сеточное разбиение для h -----------------------------------

    vector<double>x = genRegGrd(a, b, n);

    // формируем табличную функцию в узлах сетки
    // выведем соответсвующие значения
    vector<double>y(x.size());

    cout << "x" << "        " << "y" << "\n";

    cout << fixed << setprecision(6);

    for (int i = 0; i < x.size(); i++) {
        y[i] = sin(x[i]);
        cout << x[i] << " " << y[i] << "\n";
    }

    // сформируем сплайн

    Cubic_Int_Spl spl;

    // формируем его коэффициенты по исходным данным
    spl.update_spl(x, y);

    // формируем точки, которые не совпадают с узловыми
    vector<double>x_rd = spl.genRdGrd_reg();

    vector<vector<double> > results_h;


    cout << "S        " << "S'       " << "S''    " << "\n";


    // вычисляем значение сплайна и его двух первых производных в точках, которые не совпадают с узловыми
    for(int i = 0; i < x_rd.size(); i++) {

        vector<double>res = spl.get_spl_and_two_dev(x_rd[i]);

        cout << fixed << setprecision(6);

        cout << res[0] << " " << res[1] << " " << res[2] << "\n";

        results_h.push_back(res);

    }

    vector<double>err_h = spl.err_of_aprx();

    cout << "Точность сплайн-апроксимации функции sin(x) (sin(x),S;cos(x),S';-sin(x),S''): " << err_h[0] << " " <<
        err_h[1] << " " << err_h[2] << "\n";

    // формируем равномерное сеточное разбиение для h/2 -----------------------------------

    vector<double>x_h_2 = genRegGrd(a, b, n * 2);

    // формируем табличную функцию в узлах сетки
    // выведем соответсвующие значения
    vector<double>y_h_2(x_h_2.size());

    cout << "x_h_2" << "    " << "y_h_2" << "\n";

    cout << fixed << setprecision(6);

    for (int i = 0; i < x_h_2.size(); i++) {
        y_h_2[i] = sin(x_h_2[i]);
        cout << x_h_2[i] << " " << y_h_2[i] << "\n";
    }

    // формируем его коэффициенты по исходным данным

    spl.update_spl(x_h_2, y_h_2);

    // формируем точки, которые не совпадают с узловыми
    vector<double>x_rd_h_2 = spl.genRdGrd_reg();

    vector<vector<double> > results_h_2;


    cout << "S        " << "S'       " << "S''    " << "\n";


    // вычисляем значение сплайна и его двух первых производных в точках, которые не совпадают с узловыми
    for(int i = 0; i < x_rd_h_2.size(); i++) {
        vector<double>res = spl.get_spl_and_two_dev(x_rd_h_2[i]);

        cout << fixed << setprecision(6);

        cout << res[0] << " " << res[1] << " " << res[2] << "\n";

        results_h_2.push_back(res);

    }

    vector<double>err_h_2 = spl.err_of_aprx();

    cout << "Точность сплайн-апроксимации функции sin(x) (sin(x),S;cos(x),S';-sin(x),S''): " << err_h_2[0] << " " <<
        err_h_2[1] << " " << err_h_2[2] << "\n";

    // формируем равномерное сеточное разбиение для h/4 -----------------------------------

    vector<double>x_h_4 = genRegGrd(a, b, n * 4);

    // формируем табличную функцию в узлах сетки
    // выведем соответсвующие значения
    vector<double>y_h_4(x_h_4.size());

    cout << "x_h_4" << "    " << "y_h_4" << "\n";

    cout << fixed << setprecision(6);

    for (int i = 0; i < x_h_4.size(); i++) {
        y_h_4[i] = sin(x_h_4[i]);
        cout << x_h_4[i] << " " << y_h_4[i] << "\n";
    }

    // формируем его коэффициенты по исходным данным

    spl.update_spl(x_h_4, y_h_4);

    // формируем точки, которые не совпадают с узловыми
    vector<double>x_rd_h_4 = spl.genRdGrd_reg();

    vector<vector<double> > results_h_4;


    cout << "S        " << "S'       " << "S''    " << "\n";


    // вычисляем значение сплайна и его двух первых производных в точках, которые не совпадают с узловыми
    for(int i = 0; i < x_rd_h_4.size(); i++) {
        vector<double>res = spl.get_spl_and_two_dev(x_rd_h_4[i]);

        cout << fixed << setprecision(6);

        cout << res[0] << " " << res[1] << " " << res[2] << "\n";

        results_h_4.push_back(res);

    }

    vector<double>err_h_4 = spl.err_of_aprx();

    cout << "Точность сплайн-апроксимации функции sin(x) (sin(x),S;cos(x),S';-sin(x),S''): " << err_h_4[0] << " " <<
        err_h_4[1] << " " << err_h_4[2] << "\n";

    // ПУНКТ 4

    // рассмотрим наше исходное разбиение x

    spl.update_spl(x, y);

    vector<vector<double> >res_for_spl;
    vector<vector<double> >res_for_fd;
    vector<double>res_l;

    for(int i = 0; i < x.size() - 2; i++) {
        res_l = spl.get_spl_and_two_dev(x[i]);
        res_for_spl.push_back(res_l);

        // реализуем вычисление производной посредством метода конечных разностей
        // формула для вычисления первой производной с первым порядком точности
        double f_d = (sin(x[i + 1]) - sin(x[i])) / (x[1] - x[0]);
        // формула для вычисления второй производной с первым порядком точности
        double f_2d = (sin(x[i]) - 2.0 * sin(x[i + 1]) + sin(x[i + 2])) / pow(x[1] - x[0], 2);

        vector<double>f_data = {f_d, f_2d};

        res_for_fd.push_back(f_data);

    }

    // правосторонняя с первым порядком точности
    double f_d_n_1 = (sin(x[x.size() - 2]) - sin(x[x.size() - 3])) / (x[1] - x[0]);
    double f_d_n_2 = (sin(x[x.size() - 1]) - sin(x[x.size() - 2])) / (x[1] - x[0]);


    double f_d2_n_1 = (sin(x[x.size() - 2]) - 2.0 * sin(x[x.size() - 3]) + sin(x[x.size() - 4])) / pow(x[1] - x[0], 2);
    double f_d2_n_2 = (sin(x[x.size() - 1]) - 2.0 * sin(x[x.size() - 2]) + sin(x[x.size() - 3])) / pow(x[1] - x[0], 2);

    res_for_fd.push_back({f_d_n_1, f_d2_n_1});

    res_for_fd.push_back({f_d_n_2, f_d2_n_2});

    res_l = spl.get_spl_and_two_dev(x[x.size() - 2]);
    res_for_spl.push_back(res_l);
    res_l = spl.get_spl_and_two_dev(x[x.size() - 1]);
    res_for_spl.push_back(res_l);

    cout << "f'        " << "f''      " << "S'       " << "S''     " << "\n";

    for(int i = 0; i < res_for_fd.size(); i++) {

        cout << fixed << setprecision(6);

        cout << res_for_fd[i][0] << " " << res_for_fd[i][1] << " " << res_for_spl[i][1] << " " << res_for_spl[i][2] <<
            "\n";

    }


    vector<double>err_spl = spl.err_of_aprx();

    cout << "Точность сплайн-апроксимации функции sin(x) (sin(x),S;cos(x),S';-sin(x),S''): " << err_spl[0] << " " <<
        err_spl[1] << " " << err_spl[2] << "\n";


    vector<vector<double> >res_for_fd_err;

    // в рамках данного ПЗ будем вычислять погрешность используя норму пространства C[a,b]

    int N = 100000;

    double err_d;
    double err_d2;
    double max_err_fd_dev = 0;
    double max_err_fd_dev2 = 0;

    vector<double>grd = genRegGrd(a, b, N);

    for(int i = 0; i < N - 2; i++) {

        // формула для вычисления первой производной с первым порядком точности
        double f_d = (sin(grd[i + 1]) - sin(grd[i])) / (grd[1] - grd[0]);
        // формула для вычисления второй производной с первым порядком точности
        double f_2d = (sin(grd[i]) - 2.0 * sin(grd[i + 1]) + sin(grd[i + 2])) / pow(grd[1] - grd[0], 2);

        err_d = fabs(f_d - cos(grd[i]));
        err_d2 = fabs(f_2d - (-sin(grd[i])));

        if(max_err_fd_dev < err_d) {
            max_err_fd_dev = err_d;
        }
        if(max_err_fd_dev2 < err_d2) {
            max_err_fd_dev2 = err_d2;
        }

    }

    // правосторонняя с первым порядком точности
    double f_d_N_1 = (sin(grd[grd.size() - 2]) - sin(grd[grd.size() - 3])) / (grd[1] - grd[0]);
    double f_d_N_2 = (sin(grd[grd.size() - 1]) - sin(grd[grd.size() - 2])) / (grd[1] - grd[0]);

    err_d = fabs(f_d_N_1 - cos(grd[grd.size() - 2]));
    if(max_err_fd_dev < err_d) {
        max_err_fd_dev = err_d;
    }
    err_d = fabs(f_d_N_2 - cos(grd[grd.size() - 1]));
    if(max_err_fd_dev < err_d) {
        max_err_fd_dev = err_d;
    }

    double f_d2_N_1 = (sin(grd[grd.size() - 2]) - 2.0 * sin(grd[grd.size() - 3]) + sin(grd[grd.size() - 4])) /
                      pow(grd[1] - grd[0], 2);
    double f_d2_N_2 = (sin(grd[grd.size() - 1]) - 2.0 * sin(grd[grd.size() - 2]) + sin(grd[grd.size() - 3])) /
                      pow(grd[1] - grd[0], 2);

    err_d2 = fabs(f_d2_N_1 - (-sin(grd[grd.size() - 2])));
    if(max_err_fd_dev2 < err_d2) {
        max_err_fd_dev2 = err_d2;
    }
    err_d2 = fabs(f_d2_N_2 - (-sin(grd[grd.size() - 1])));
    if(max_err_fd_dev2 < err_d2) {
        max_err_fd_dev2 = err_d2;
    }

    vector<double>err_fd = {
        max_err_fd_dev, max_err_fd_dev2
    };

    cout << "Точность аппроксимации с помощью МКР функции sin(x) (cos(x),f';-sin(x),f''): " << err_fd[0] << " " <<
        err_fd[1] << "\n";



    return 0;
}
