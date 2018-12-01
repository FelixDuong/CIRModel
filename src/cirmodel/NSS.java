/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cirmodel;

import java.io.IOException;
import java.util.Arrays;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import static org.apache.commons.math3.util.MathArrays.sequence;

import org.ejml.simple.SimpleMatrix;


import java.awt.Color;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;


/**
 *
 * @author tuanda
 */
public class NSS {

    static class LineChartDemo2 extends ApplicationFrame {
/**
* Creates a new demo.
*
* @param title the frame title.
*/
LineChartDemo2(String title) {
super(title);
XYDataset dataset = createDataset();
JFreeChart chart = createChart(dataset);
ChartPanel chartPanel = new ChartPanel(chart);
chartPanel.setPreferredSize(new java.awt.Dimension(700, 500));
setContentPane(chartPanel);
}
/**
* Creates a sample dataset.
*
* @return a sample dataset.
*/
private static XYDataset createDataset() {
XYSeries series1 = new XYSeries("Yields");
for(int i=0; i<Y.getNumElements();i++){
    series1.add(i+1,Y.get(i, 0));
}
XYSeriesCollection dataset = new XYSeriesCollection();
dataset.addSeries(series1);

return dataset;
}
/**
* Creates a chart.
*
* @param dataset the data for the chart.
*
* @return a chart.
*/
private static JFreeChart createChart(XYDataset dataset) {
// create the chart...
JFreeChart chart = ChartFactory.createXYLineChart(
"Yield Curve", // chart title
"Months", // x axis label
"%", // y axis label
dataset, // data
PlotOrientation.VERTICAL,
true, // include legend
true, // tooltips
false // urls
);
// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
chart.setBackgroundPaint(Color.white);
// get a reference to the plot for further customisation...
XYPlot plot = (XYPlot) chart.getPlot();
plot.setBackgroundPaint(Color.lightGray);
plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
plot.setDomainGridlinePaint(Color.white);
plot.setRangeGridlinePaint(Color.white);
XYLineAndShapeRenderer renderer
= (XYLineAndShapeRenderer) plot.getRenderer();
renderer.setShapesVisible(true);
renderer.setShapesFilled(true);
// change the auto tick unit selection to integer units only...
NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
// OPTIONAL CUSTOMISATION COMPLETED.
return chart;
}
/**
* Creates a panel for the demo (used by SuperDemo.java).
*
* @return A panel.
*/
private static JPanel createDemoPanel() {
JFreeChart chart = createChart(createDataset());
return new ChartPanel(chart);
}
/**
* Starting point for the demonstration application.
*
* @param args ignored.
*/

}
    
    public static SimpleMatrix T, Y;

    private static SimpleMatrix NSS(double[] x, SimpleMatrix T) {

        SimpleMatrix one = new SimpleMatrix(T.getNumElements(), 1);
        one.fill(1.0);
        double lambda_1 = 1. / x[4];
        double lambda_2 = 1. / x[5];
        SimpleMatrix f, T_lambda_1, T_lambda_2, T_lambda_1_exp, T_lambda_2_exp, temp1, temp2, one_1, one_2;
        T_lambda_1 = T.scale(lambda_1);
        T_lambda_2 = T.scale(lambda_2);
        T_lambda_1_exp = T.scale(-1 * lambda_1).elementExp();
        T_lambda_2_exp = T.scale(-1 * lambda_2).elementExp();

        one_1 = one.minus(T_lambda_1_exp);
        one_2 = one.minus(T_lambda_2_exp);
        temp1 = one_1.elementDiv(T_lambda_1);
        f = temp1.scale(x[1]).plus(x[0]);
        temp2 = one_2.elementDiv(T_lambda_2);

        temp1 = temp1.minus(T_lambda_1_exp).scale(x[2]);
        temp2 = temp2.minus(T_lambda_2_exp).scale(x[3]);

        // f = T.divide(-1*x[4]).elementExp()
        f = f.plus(temp1).plus(temp2);

        return f;
    }

    private static double NSS_func(double[] x) {

        SimpleMatrix one = new SimpleMatrix(10, 1, true, new double[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
        double lambda_1 = 1. / x[4];
        double lambda_2 = 1. / x[5];
        SimpleMatrix f, T_lambda_1, T_lambda_2, T_lambda_1_exp, T_lambda_2_exp, temp1, temp2, one_1, one_2;
        T_lambda_1 = T.scale(lambda_1);
        T_lambda_2 = T.scale(lambda_2);
        T_lambda_1_exp = T.scale(-1 * lambda_1).elementExp();
        T_lambda_2_exp = T.scale(-1 * lambda_2).elementExp();

        one_1 = one.minus(T_lambda_1_exp);
        one_2 = one.minus(T_lambda_2_exp);
        temp1 = one_1.elementDiv(T_lambda_1);
        f = temp1.scale(x[1]).plus(x[0]);
        temp2 = one_2.elementDiv(T_lambda_2);

        temp1 = temp1.minus(T_lambda_1_exp).scale(x[2]);
        temp2 = temp2.minus(T_lambda_2_exp).scale(x[3]);

        // f = T.divide(-1*x[4]).elementExp()
        f = f.plus(temp1).plus(temp2).minus(Y).elementPower(2);

        return f.elementSum();
    }

    public static void main(String args[]) {

        //double[] x = new double[]{4.143051872,8.12317E-05,197591.2447,-1156978.376,5988.928358,35301.06808};
        T = new SimpleMatrix(10, 1, true, new double[]{0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20});
        Y = new SimpleMatrix(10, 1, true, new double[]{3.75, 3.78, 3.85, 3.99, 4.13, 4.39, 4.62, 4.95, 5.40, 5.74});
        try{
            //Y.saveToFileCSV("input_yields.csv");
            Y = Y.loadCSV("input_yields.csv");
        //System.out.println("Loss function NSS :  " + NSS_func(x));
        }
        catch (IOException e) {
        throw new RuntimeException(e);
         }
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-12, 1e-15);
        //PowellOptimizer optimizer = new PowellOptimizer(1e-9, 1e-15);
        NSS_objfunc objfunc = new NSS_objfunc();
        PointValuePair optimum;
        optimum = optimizer.optimize(
                new MaxEval((int) 1e12),
                new MaxIter(5000),
                new ObjectiveFunction(objfunc),
                GoalType.MINIMIZE,
                new InitialGuess(new double[]{5, -2, -1, -2, 2, 3}),
                new NelderMeadSimplex(6));
        // new SimpleBounds(new double[]{0.0001,-1000,-1000,-1000,0.0001,0.0001},new double[]{1000,1000,1000,1000,1000,1000}));
//                          //new MultiDirectionalSimplex(new double[]{alpha, mu, sigma}));
//
        System.out.println("Optimization point: " + Arrays.toString(optimum.getPoint()));
        System.out.println("Loss function NSS optimized :  " + NSS_func(optimum.getPoint()));

       
        
        
        // Generate yield series from the result
        int[] t = sequence(240, 1, 1);

        SimpleMatrix temp = new SimpleMatrix(t.length, 1);

        //System.out.println("Yield curve : " +  NSS(optimum.getPoint(),temp).toString());
        for (int i = 0; i < t.length; ++i) {
            temp.set(i, 0, t[i] * 0.083333);
        }
        try {
            temp = NSS(optimum.getPoint(), temp);
            temp.reshape(t.length, 1);
            temp.saveToFileCSV("NSS_Yields.csv");
            Y = temp;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        
        LineChartDemo2 demo = new LineChartDemo2("Yield Curve");
        demo.pack();
        RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);

    }

    private static class NSS_objfunc implements MultivariateFunction {

        @Override
        public double value(double[] point) {
            return NSS_func(point);
        }

    }

}
