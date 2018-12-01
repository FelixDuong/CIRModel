/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cirmodel;

/**
 *
 * @author felixduong
 */

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author felixduong
 */
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import java.io.IOException;
import java.io.Reader;
import static java.lang.Double.isNaN;
import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import static org.apache.commons.math3.util.MathArrays.ebeDivide;
import static org.apache.commons.math3.util.MathArrays.ebeSubtract;
import umontreal.ssj.probdist.ChiSquareNoncentralDist;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import static org.apache.commons.math3.stat.StatUtils.min;
import static org.apache.commons.math3.stat.StatUtils.sum;
import static org.apache.commons.math3.util.MathArrays.checkNotNaN;
import static org.apache.commons.math3.util.MathArrays.distance1;
import static org.apache.commons.math3.util.MathArrays.scale;
import org.ejml.simple.SimpleMatrix;

public class CIRModel {

    private static final String SAMPLE_CSV_FILE_PATH = "data.csv";
    public static double[] r_, r_0, r_1;

//    private static double loglikeCIR(double[] params) { // Vie RemiRelard
//        double alpha, beta, sigma, q1, q2;
//        double h = 1. / 360.;
//        SimpleMatrix tau = new SimpleMatrix(r_.length, 1);
//        SimpleMatrix one = new SimpleMatrix(r_.length, 1);
//        tau.fill(3.0);// tenor selected
//        one.fill(1.0);
//        alpha = params[0];
//        beta = params[1] / 100;
//        sigma = params[2] / 10;
////        q1 = params[3] / 10;
////        q2 = params[4] * 10;
//        q2 = 10;
//        q1 = 0;
//        double a = alpha + q2 * sigma;
//
//        double b = (alpha * beta - q1 * sigma) / a;
//
//        double gam = sqrt(a * a + 2 * sigma * sigma);
//
//        SimpleMatrix d = one.minus(tau.scale(-1 * gam).elementExp());
//
//        SimpleMatrix B = d.elementDiv(d.scale(0.5 * (gam + a)).plus(one.minus(d).scale(gam)));
//
//        SimpleMatrix A = (tau.scale(0.5 * (gam - a)).plus((d.scale(0.5 * (1 + a / gam)).plus(one.minus(d))).elementLog())).scale(-(2 * a * b / (sigma * sigma)));
//
//        SimpleMatrix LL = new SimpleMatrix(r_.length - 1, 1);
//        LL.fill(1.0);
//
//        if (a <= 0 || b <= 0) {
//            LL = LL.scale(1E20);
//            return (log(LL.elementSum()));
//        }
//        SimpleMatrix r_01 = new SimpleMatrix(r_.length, 1);
//        for (int i = 0; i < r_.length; ++i) {
//            r_01.set(i, 0, r_[i]);
//        }
//        r_01 = tau.elementMult(r_01).plus(A.scale(100)).elementDiv(B); // scale r serie orginal 
//        for (int i = 0; i < r_.length; ++i) {
//            if (r_01.get(i, 0) <= 0) {
//                LL = LL.scale(1E20);
//                return (log(LL.elementSum()));
//            }
//        }
//
//        //parameters for the returns
//        double phi = exp(-alpha * h);
//        double nu = 4 * alpha * beta / (sigma * sigma);
//        double omega = beta * (1 - phi) / nu;
//        r_01 = r_01.scale(1 / omega);
//        SimpleMatrix D = r_01.scale(phi);
//        D = D.extractMatrix(1, D.getNumElements(), 0, 1);
//        B = B.extractMatrix(1, B.getNumElements(), 0, 1);
//        for (int i = 0; i < B.getNumElements(); ++i) {
//            B.set(i, 0, abs(B.get(i, 0)));
//        }
//        tau = tau.extractMatrix(1, tau.getNumElements(), 0, 1);
//        r_01 = r_01.extractMatrix(1, r_01.getNumElements(), 0, 1);
//
//        ChiSquareNoncentralDist xi;
//
//        SimpleMatrix gpdf = new SimpleMatrix(r_01.getNumElements(), 1);
//        //gpdf.zero();
//        for (int t = 0; t < r_01.getNumElements(); t++) {
//
//            xi = new ChiSquareNoncentralDist(nu, D.get(t, 0));
//            double den_ = max(xi.density(r_01.get(t, 0)), 1E-15);
//            gpdf.set(t, 0, den_);
//
//            //gpdf.addToEntry(i, max(xi.density(s.getEntry(i)), eps));
//            //gpdf.addToEntry(i, max(xi.cdf(s.getEntry(i)),eps));
//        }
//
//        LL = B.elementDiv(tau).elementLog().plus(log(omega)).minus(gpdf.elementLog());
//
//        for (int i = 0; i < LL.getNumElements(); ++i) {
//            if (LL.get(i, 0) <= 0 || isNaN(LL.get(i, 0))) {
//                LL = LL.scale(1E20);
//                return log((LL.elementSum()));
//            }
//
//        }
//        return log((LL.elementSum()));
//    }
    private static double loglikeCIR(double[] params) {

        double alpha, mu, sigma;
        alpha = params[0];
        mu = params[1];
        sigma = params[2];
        //SimpleMatrix DataF, DataL, gpdf;
        double[] DataF, DataL, gpdf, u, v, s, nc;
        double TimeStep = 1. / 360.;
//            DataF = new SimpleMatrix(r_1.length, 1, true, r_1);
//            DataL = new SimpleMatrix(r_0.length, 1, true, r_0);
        DataF = r_1;
        DataL = r_0;
        u = new double[r_0.length];
        v = new double[r_0.length];
        s = new double[r_0.length];
        nc = new double[r_0.length];
        gpdf = new double[r_0.length];
        double c = 2 * alpha / (sigma * sigma * (1 - exp(-alpha * TimeStep)));
        double q = 2 * alpha * mu / (sigma * sigma) - 1;
//            SimpleMatrix u = DataL.scale(c * exp(-alpha * TimeStep));
//            SimpleMatrix v = DataF.scale(c);
//            SimpleMatrix s = DataF.scale(2 * c);
//            SimpleMatrix nc = u.scale(2); // noncentrality parameter
        for (int i = 0; i < r_0.length; i++) {
            u[i] = c * exp(-alpha * TimeStep) * DataL[i];
            s[i] = 2 * c * DataF[i];
            nc[i] = 2 * u[i];

        }
        double df = 2 * q + 2;
        ChiSquareNoncentralDist xi;
        //double[] gpdf= new double[r_0.length];
        // gpdf = new double[r_0.length];
        // gpdf = new SimpleMatrix(r_1.length, 1);
        for (int i = 0; i < r_0.length; ++i) {
            //if (df > 0 && nc.getEntry(i) > 0) {
            xi = new ChiSquareNoncentralDist(df, nc[i]);
            //gpdf.addToEntry(i, xi.density(s.getEntry(i)));
            // System.out.println("\n density func value : " + xi.density(s.getEntry(i)) + "\n");
            //System.out.println("\n density value: " + xi.density(s[i]));
            gpdf[i] = -log(2 * c * xi.density(s[i]));
            //gpdf.addToEntry(i, max(xi.cdf(s.getEntry(i)),eps));
            // } else {
            //     gpdf.addToEntry(i, 0);
            // }

        }
        double LL = sum(gpdf);
        System.out.println("\nlog-likelihood value : " + LL + "\n");
        return LL;
    }

    private static double[] gradient_log(double[] params) {
        final double eps = 0.0001;

        double[] grd, x1, x2;
        x1 = params.clone();
        x2 = params.clone();
        grd = new double[]{0, 0, 0};
        for (int i = 0; i < params.length; ++i) {
            x1[i] += eps;
            x2[i] -= eps;
            // if (x1[i] == 0 || x2[i] == 0) {
            //    grd[i] = 0;
            //   } else {
            grd[i] = (loglikeCIR(x1) - loglikeCIR(x2)) / (2 * eps);
            // }
        }
        return grd;

    }

    private static double[] SecantMethod(double[] params) {
        double eps = 1E-10;
        double fn_1, fn;
        double[] xn_1, xn, temp;
        double[] Q, Q_0, delta_Q;
        double ratio = 0.001;
        xn_1 = params;
        xn = new double[]{params[0] * 0.9, params[1] * 0.9, params[2] * 0.9};
        Q = new double[]{0, 0, 0};
        System.out.println("Optimization start: " + loglikeCIR(xn_1));
        double delta_f;
        fn_1 = loglikeCIR(xn_1);
        fn = loglikeCIR(xn);
        delta_f = fn_1 - fn;
        // Q = new double[]{delta_f / (params[0] * 0.5), delta_f / (params[1] * 0.5), delta_f / (params[2] * 0.5)};
        // Q = new double[]{0.01, 0.01, 0.01};
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < params.length; j++) {
                xn[j] = xn_1[j] - Q[j];
                fn_1 = loglikeCIR(xn_1);
                fn = loglikeCIR(xn);
                if ((fn_1 < fn) || xn[j] > xn_1[j]) {
                    //xn[j] = xn_1[j] + Q[j];
                    Q[j] = delta_f / (-xn_1[j] + xn[j]);
                    fn_1 = loglikeCIR(xn_1);
                    fn = loglikeCIR(xn);
                }
                delta_f = fn_1 - fn;
                if (delta_f < eps && delta_f > 0) {
                    System.out.println("Optimization value: " + loglikeCIR(xn));
                    return xn;
                }
                Q[j] = delta_f / (xn_1[j] - xn[j]);
                if (xn[j] > 0) {
                    xn_1[j] = xn[j];
                }
                if (isNaN(Q[j])) {
                    Q[j] = 0.001;
                }

            }
//            System.out.println("Optimization value: " + loglikeCIR(xn));
//            return xn;
        }
        System.out.println("Optimization value: " + loglikeCIR(xn));
        return xn;
    }

    private static double[] NewtonRaphson(double[] params) {
        double eps = 0.0000001;
        double dis, dis_h, log_target;
        System.out.println("prams pass in : " + Arrays.toString(params));
        System.out.println("loglikehood : " + loglikeCIR(params));
        System.out.println("Gradient of params : " + Arrays.toString(gradient_log(params)));
        log_target = loglikeCIR(params);

        checkNotNaN(gradient_log(params));
        double[] update_params, h_, h;
        h = new double[]{loglikeCIR(params), loglikeCIR(params), loglikeCIR(params)};
        //h = ebeDivide(h, gradient_log(params));
        h = scale(-1, ebeDivide(h, gradient_log(params)));
        h_ = h.clone();
        System.out.println("h : " + Arrays.toString(h));
        update_params = params.clone();
        dis = log_target - loglikeCIR(params);
        //run search
        for (int i = 0; i < 100; ++i) {
            for (int j = 0; j < 3; j++) {
                params[j] += h[j];
                checkNotNaN(params);
                if (loglikeCIR(params) < log_target && params[j] > 0 && params[0] < 10 && !Double.isInfinite(loglikeCIR(params))) {
                    update_params[j] = params[j];
                    dis = log_target - loglikeCIR(params);
                    log_target = loglikeCIR(params);
                } else {
                    System.out.println("Reverse direction params no : " + j);
                    h[j] = -1 * h[j];
                    params[j] = update_params[j];
                }

            }

            System.out.println("=========================");
            System.out.println("Iterations No: " + i);
            System.out.println("Params changes: " + Arrays.toString(update_params));

            System.out.println("Gradient of params changes: " + Arrays.toString(gradient_log(params)));
            System.out.println("Objective value loglike: " + loglikeCIR(params));

            //dis = log_target - loglikeCIR(params);
            params = update_params.clone();
            System.out.println("IS INFINITY LOG FUNC : " + Double.isInfinite(loglikeCIR(params)));
            checkNotNaN(gradient_log(params));
            //h = ebeDivide(h, gradient_log(params));
            h = new double[]{loglikeCIR(params), loglikeCIR(params), loglikeCIR(params)};
            h = ebeDivide(h, gradient_log(params));
            dis_h = distance1(h, h_);
            checkNotNaN(h);

            System.out.println("h new change: " + Arrays.toString(h));
            System.out.println("Distance new h: " + dis_h);
            System.out.println("Distance loglikehood: " + dis);
            System.out.println("=========================");
//            if (dis < eps ) {
//                h = new double[]{loglikeCIR(params), loglikeCIR(params), loglikeCIR(params)};
//                h = ebeDivide(h, gradient_log(params));
//            }
//
            if (abs(dis) < eps && i > 10) {
                return params;
            }

            h_ = h.clone();
            log_target = loglikeCIR(params);

        }
        return params;
    }

    private static double[] NelderMead(double[] params, int n) {
        double alpha0 = 1;
        double beta0 = 1;
        double gamma0 = 0.5;
        double delta0 = 0.5;

        double[] f = new double[n + 1];
        double[][] x = new double[n + 1][params.length];
        double[] params_r = new double[params.length]; // params to pass loglikeCIR to cacul f value
        double[] params_n_1 = new double[params.length];
        //double[] x = new double[n+1]; // serie x of parameter
        //double[] f = new double[n+1]; // value function f stored here
        for (int j = 0; j <= n; ++j) {
            x[0] = params.clone();
            for (int i = 0; i < params.length; ++i) {
                if (j > 0) {
                    x[j][i] = x[j - 1][i] - 0.0025;
                }
            }
            f[j] = loglikeCIR(x[j]);
        }
        // Sort process series x[]
        int m = n;
        do {
            for (int i = n - m; i < n; i++) {
                if (f[i] == min(f, n - m, m)) {
                    double[] temp = x[i];
                    x[i] = x[n - m];
                    x[n - m] = temp;
                    double f_temp = f[i];
                    f[i] = f[n - m];
                    f[n - m] = f_temp;

                }
            }
            m -= 1;

        } while (m == 0);

        System.out.println("x serie : " + Arrays.deepToString(x));
        System.out.println("f serie : " + Arrays.toString(f));
        System.out.println("Loglike value after adjusted : " + loglikeCIR(x[0]));
//            
//            double x_mean = mean(copyOfRange(x,0,n-1));
//            double x_r = x_mean + alpha0*(x_mean - x[n]);
//            params_r[i]= x_r;
//            params_n_1[i] = x[n];
//            if(loglikeCIR(params_r)>= loglikeCIR(params_0) && loglikeCIR(params_r)< loglikeCIR(params_r) ){
//                x[n] = x_r;

        // Chua sort chuoi fx , va order x
//        }
        return params;
    }

    @SuppressWarnings("empty-statement")
    public static void main(String[] args) throws IOException {
        int i = 0;
        try (
                Reader reader = Files.newBufferedReader(Paths.get(SAMPLE_CSV_FILE_PATH));
                CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT);) /* // Reading all records at once into memory
            List<CSVRecord> csvRecords = csvParser.getRecords();
            But you should avoid this method if you’re reading a significantly large CSV file. 
            You might run into memory issues because the getRecords() method loads the entire CSV contents into memory.*/ {
            double h = 0.00278;
            double[] diff;
            double[][] regressors;
            RealMatrix A;
            RealVector B;
            r_ = new double[5000];
            //r_0 = new double[1107];
            //r_1 = new double[1107];
            //diff = new double[1107];
            for (CSVRecord csvRecord : csvParser) {
                // Accessing Values by Column Index
                String temp = csvRecord.get(0);
                r_[i] = Double.parseDouble(temp);
                //System.out.println(r_[i]);
                ++i;

            }
            r_ = Arrays.copyOf(r_, i);
            //System.out.println(mean(r_));
            r_0 = Arrays.copyOfRange(r_, 0, r_.length - 1);
            r_1 = Arrays.copyOfRange(r_, 1, r_.length);
            diff = ebeSubtract(r_0, r_1);
            regressors = new double[r_0.length][2];

            for (i = 0; i < r_0.length; ++i) {
                r_0[i] = pow(r_0[i], 0.5);
                regressors[i][0] = h / r_0[i];
                regressors[i][1] = h * r_0[i];

            }
            diff = ebeDivide(diff, r_0);
            A = new Array2DRowRealMatrix(regressors, false);
            B = new ArrayRealVector(diff, false);

            DecompositionSolver solver = new SingularValueDecomposition(A).getSolver();
            RealVector solution = solver.solve(B);
            Array2DRowRealMatrix drift = new Array2DRowRealMatrix(new double[][]{{solution.getEntry(0), solution.getEntry(1)}}, false);

            RealMatrix res = A.multiply(drift.transpose());
            //B = B.mapMultiply(-1);

            B = res.getColumnVector(0).subtract(B);

            double alpha = drift.getEntry(0, 1);
            double mu = -drift.getEntry(0, 0) / drift.getEntry(0, 1);
            double sigma = new Variance().evaluate(B.toArray());
            sigma = sqrt(sigma / h);
            r_0 = Arrays.copyOfRange(r_, 0, r_.length - 1);
            System.out.println("Alpha :" + alpha);
            System.out.println("mu - Beta : " + mu);
            System.out.println("Sigma :" + sigma);
            System.out.println("Loglikelihood : " + loglikeCIR(new double[]{alpha, mu, sigma}));

//System.out.println("Loglikelihood test : " + loglikeCIR(new double[]{6430.760175368759, 3.426435379233, 0.5955618328539553}));
            // System.out.println("Gradient of Loglikelihood : " + Arrays.toString(gradient_log(new double[]{alpha, mu, sigma})));
            //Optimization
//            double[] rls = NelderMead(new double[]{alpha, mu, sigma},10);
//            double[] rls = NewtonRaphson(new double[]{alpha, mu, sigma});
//            double[] rls = SecantMethod(new double[]{2, 1, 1});
//            System.out.println("Optimization parameters: " + Arrays.toString(rls));
//            BOBYQAOptimizer optimizer = new BOBYQAOptimizer(7, 100, 1E-12);
//            logfunc objfunc = new logfunc();
//            PointValuePair optimum
//                    = optimizer.optimize(
//                            new MaxEval((int) 1e10),
//                            new ObjectiveFunction(objfunc),
//                            GoalType.MINIMIZE,
//                            new InitialGuess(new double[]{alpha, mu, sigma,0,0}),
//                            //new SimpleBounds(boundaries[0], boundaries[1]));
//                            new SimpleBounds(new double[]{0.0, 0.0, 0.1,0,0},
//                            new double[]{100, 100, 100,100,100}));
            SimplexOptimizer optimizer = new SimplexOptimizer(1e-5, 1e-10);
            // PowellOptimizer optimizer = new PowellOptimizer(1e-3, 1e-10);
            logfunc objfunc = new logfunc();
            PointValuePair optimum
                    = optimizer.optimize(
                            new MaxEval((int) 1E10),
                            new MaxIter(1000),
                            new ObjectiveFunction(objfunc),
                            GoalType.MINIMIZE,
                            new InitialGuess(new double[]{alpha, mu, sigma}),
                            new NelderMeadSimplex(3));
//                            new MultiDirectionalSimplex(new double[]{alpha, mu, sigma}));
//
            System.out.println(Arrays.toString(optimum.getPoint()) + " : " + Arrays.toString(optimum.getPointRef()));
            if (optimum.getPointRef()[2] < 0) {
                optimum.getPointRef()[2] = sigma;
            }
            System.out.println("Adjusted params : " + Arrays.toString(optimum.getPoint()));
            System.out.println("Objective value : " + loglikeCIR(optimum.getPoint()));
////             //Simulation
            Sim_CIR run = new Sim_CIR();
//            run.run_sim(rls, 4.1367, 100, 90);
            run.run_sim(optimum.getPoint(), 4.1367, 100, 90);
//
        } catch (IOException ex) {
            System.out.println(ex.toString());
        }
    }

    private static class logfunc implements MultivariateFunction {

        @Override
        public double value(double[] point) {
            return loglikeCIR(point);
        }

    }

//    private static class NSS implements MultivariateFunction {
//        private static double[] T;
//        @Override
//        public double value(double[] point) {
//            
//            double y = point[0] + point[1] * (1 - exp(-T / point[4])) * 1 / (T / point[4])
//			+ point[2] * ((1 - exp(-T / point[4])) * 1 / (T / point[4]) - exp(-T / point[4]))
//			+ point[3] * ((1 - exp(-T / point[5])) * 1 / (T / point[5]) - exp(-T / point[5]));
//	return y;
//            
//            
//            
//        }
//
//    }
}
