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
import org.apache.commons.csv.CSVPrinter;

import java.io.*;
import static java.lang.Math.exp;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;
import umontreal.ssj.randvar.ChiSquareNoncentralGamGen;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.rng.LFSR113;

public class Sim_CIR {
    private static final String SAMPLE_CSV_FILE = "CIRmodel_sim.csv";
    RandomVariateGen gen;
    public void run_sim(double[] params, double r0, int sim_no, int length_time) {
        
        double sigma_, alpha_, beta_,h;
        sigma_ = params[2];
        alpha_ = params[0];
        beta_ = params[1];
        h = 1./360.;
        SimpleMatrix r = new SimpleMatrix(length_time,sim_no);
           
        
//        sigmaSquared = sigma^2;
//nu           = 4 * alpha * beta / sigmaSquared 
//phi          = exp( - alpha * h );
//omega        = sigmaSquared * ( 1-phi )  / ( 4 * alpha );
//
//r    = zeros(n+1,1);
//r(1) = r0;
//
//addpath( '../ChapterA')
//for t = 2:n+1
//    x = r(t-1) / omega;
//    D = x * phi;  % non-centrality parameter
//    r(t) = omega * SimNChi2( nu, D );
//   % r(t) = omega * ncx2rnd( nu, D );
//end
                
        double	sigmaSquared = sigma_ * sigma_;
	double	nu = 4 * alpha_ * beta_ / sigmaSquared;
	double  phi = exp(-alpha_ * h);
	double  omega = sigmaSquared * (1 - phi) / (4 * alpha_);
           
        
        for (int i = 0; i < sim_no; ++i) {
		r.set(0,i,r0);
		
		for (int t = 1; t < length_time; ++t) {
			double x = r.get(t-1, i) / omega;
			double D = x * phi;
			RandomVariateGen gen = new ChiSquareNoncentralGamGen (new LFSR113(), nu, D);
			r.set(t,i,omega*generate(gen));
		}
        }
	try{	
        r.saveToFileCSV(SAMPLE_CSV_FILE);
        
//      BufferedWriter writer = null;
//            try {
//                  CSVPrinter csvPrinter;
//                writer = Files.newBufferedWriter(Paths.get(SAMPLE_CSV_FILE));
//                 csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT);
//                        
//            for (i = 0; i<sim_no; i++){
//                for (int j=0; j<length_time; j++){
//                    csvPrinter.print(r.get(j,i));
//                }
//                csvPrinter.println();
//            }
           
            
//            csvPrinter.printRecord(Arrays.asList(r));

//            csvPrinter.flush();  
            
            
            } catch (IOException ex) {
                Logger.getLogger(Sim_CIR.class.getName()).log(Level.SEVERE, null, ex);
            }

      
        
    }
    private static double generate(RandomVariateGen gen) {
      
		
	double u;
        u = gen.nextDouble();
        return u;
    
    }
}


