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
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import java.io.IOException;
import java.io.Reader;
import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import static java.util.Arrays.copyOfRange;
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
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import static org.apache.commons.math3.stat.StatUtils.mean;
import static org.apache.commons.math3.stat.StatUtils.min;
import static org.apache.commons.math3.stat.StatUtils.min;
import static org.apache.commons.math3.stat.StatUtils.min;
import static org.apache.commons.math3.util.MathArrays.checkNotNaN;
import static org.apache.commons.math3.util.MathArrays.copyOf;
import static org.apache.commons.math3.util.MathArrays.distance1;
import static org.apache.commons.math3.util.MathArrays.ebeAdd;
import static org.apache.commons.math3.util.MathArrays.scale;

public class Est_CIR {

    private static final String SAMPLE_CSV_FILE_PATH = "data.csv";
    public static double[] r_0, r_1;

    private static double loglikeCIR(double[] params) {
        final double eps = 1E-7;
        double alpha, mu, sigma;
        alpha = params[0];
        mu = params[1];
        sigma = params[2];
        ArrayRealVector DataF, DataL;
        double LL = 0;
        double TimeStep = 0.00278;
        DataF = new ArrayRealVector(r_1, false);
        DataL = new ArrayRealVector(r_0, false);
        double c = 2 * alpha / (sigma * sigma * (1 - exp(-alpha * TimeStep)));
        double q = 2 * alpha * mu / sigma * sigma - 1;
        RealVector u = DataL.mapMultiply(c * exp(-alpha * TimeStep));
        RealVector v = DataF.mapMultiply(c);
        RealVector s = DataF.mapMultiply(2 * c);
        RealVector nc = u.mapMultiply(2); // noncentrality parameter
        double df = 2 * q + 2;
        ChiSquareNoncentralDist xi;
        RealVector gpdf = s;
        for (int i = 0; i < r_0.length; ++i) {
            if (df > 0 && nc.getEntry(i) > 0) {
                xi = new ChiSquareNoncentralDist(df, nc.getEntry(i));
               gpdf.addToEntry(i, xi.density(s.getEntry(i)));
                //gpdf.addToEntry(i, max(xi.density(s.getEntry(i)), eps));
                //gpdf.addToEntry(i, max(xi.cdf(s.getEntry(i)),eps));
            } else {
                gpdf.addToEntry(i, 0);
            }

            LL += -Math.log(2 * c * gpdf.getEntry(i));
        }
        return LL;
    }

    private static double[] gradient_log(double[] params) {
        final double eps = 0.0001;

        double[] grd, x1, x2;;
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
                if (loglikeCIR(params) < log_target && params[j] > 0 && params[0]<10 && !Double.isInfinite(loglikeCIR(params))) {
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
            if (abs(dis) < eps && i>10 ) {
                return params;
            }

            h_ = h.clone();
            log_target = loglikeCIR(params);

        }
        return params;
    }
    
    private static double[] NelderMead(double[] params,int n){
        double alpha0 = 1;
        double beta0 = 1;
        double gamma0 = 0.5;
        double delta0 = 0.5;
        double temp;
        
        double[] params_r = new double[3]; // params to pass loglikeCIR to cacul f value
        double[] params_n_1 = new double[3];
        double[] x = new double[n+1]; // serie x of parameter
        double[] f = new double[n+1]; // value function f stored here
        for(int i = 0; i <params.length; ++i){
            double[] params_0 = params.clone();
            x[0] = params[i];
            for (int j =0;j<=n;++j){
                if(j>0) {x[j] = x[j-1] + 2;}
                params_0[i]= x[j];
                f[j] = loglikeCIR(params_0);
            }
            // Sort process series x[]
            for (int j=0;j<=n;++j){
                if(f[j] < min(f,j+1,n-j)) {
                temp = x[j+1]; x[j+1] = x[j]; x[j] = temp;
                temp=f[j+1]; f[j+1] = f[j]; f[j] = temp;}
            }
            System.out.println("x serie : " + Arrays.toString(x));
            System.out.println("f serie : " + Arrays.toString(f));
//            
//            double x_mean = mean(copyOfRange(x,0,n-1));
//            double x_r = x_mean + alpha0*(x_mean - x[n]);
//            params_r[i]= x_r;
//            params_n_1[i] = x[n];
//            if(loglikeCIR(params_r)>= loglikeCIR(params_0) && loglikeCIR(params_r)< loglikeCIR(params_r) ){
//                x[n] = x_r;
            }
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
            double[] r_, diff;
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
            r_0 = Arrays.copyOf(r_, r_.length - 1);
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

            System.out.println("Alpha :" + alpha);
            System.out.println("mu - Beta : " + mu);
            System.out.println("Sigma :" + sigma);
            System.out.println("Loglikelihood : " + loglikeCIR(new double[]{alpha, mu, sigma}));
            System.out.println("Loglikelihood test : " + loglikeCIR(new double[]{6430.760175368759, 3.426435379233, 0.5955618328539553}));
            // System.out.println("Gradient of Loglikelihood : " + Arrays.toString(gradient_log(new double[]{alpha, mu, sigma})));
            //Optimization
               double[] rls = NelderMead(new double[]{alpha, mu, sigma},10);
//            double[] rls = NewtonRaphson(new double[]{alpha, mu, sigma});
//            System.out.println("Optimization parameters: " + Arrays.toString(rls));
//            double x_input[] = {alpha,mu,sigma};
//            //x_input generated from Poisson(3)
//
//            MultivariateFunction f = (double[] lam) -> (loglikeCIR(x_input));
//            MultivariateOptimizer optim = new BOBYQAOptimizer(x_input.length*3);
//            PointValuePair result;
//            result = optim.optimize(new MaxEval((int) 1e10),
//                    new ObjectiveFunction(f),
//                    GoalType.MINIMIZE,
//                    new InitialGuess(new double[]{alpha,mu,sigma}) );
//            System.out.println(Arrays.toString(result.getPoint()));
//            SimplexOptimizer optimizer = new SimplexOptimizer(1e-3, 1e-10);
//            // PowellOptimizer optimizer = new PowellOptimizer(1e-3, 1e-10);
//            logfunc objfunc = new logfunc();
//            PointValuePair optimum
//                    = optimizer.optimize(
//                            new MaxEval((int) 1e10),
//                            new MaxIter(1000),
//                            new ObjectiveFunction(objfunc),
//                            GoalType.MINIMIZE,
//                            new NonNegativeConstraint(true),
//                            new InitialGuess(new double[]{alpha, mu, sigma}),
//                            new NelderMeadSimplex(3));
//            //new MultiDirectionalSimplex(new double[]{alpha, mu, sigma}));
//
//            System.out.println(Arrays.toString(optimum.getPoint()) + " : " + Arrays.toString(optimum.getPointRef()));
//            if (optimum.getPointRef()[2] < 0) {
//                optimum.getPointRef()[2] = sigma;
//            }
//            System.out.println("Adjusted params : " + Arrays.toString(optimum.getPoint()));
//            System.out.println("Objective value : " + loglikeCIR(optimum.getPoint()));
//             //Simulation
//            Sim_CIR run = new Sim_CIR();
//            run.run_sim(rls, 4.1367, 100, 90);
           // run.run_sim(optimum.getPoint(), 4.1367, 100, 90);
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

}
