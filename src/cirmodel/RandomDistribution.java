/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cirmodel;

import static java.lang.System.out;
import umontreal.ssj.probdist.ChiSquareNoncentralDist;
import umontreal.ssj.randvar.ChiSquareNoncentralGen;
import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import java.util.Arrays;
import umontreal.ssj.rng.RandomStream;
/**
 *
 * @author tuanda
 */

public class RandomDistribution {
    public double generate(RandomVariateGen gen) {
		double u;
		
		return	u = gen.nextDouble();
			
	}
//    public static void main(String[] args){
//    ChiSquareNoncentralDist xi;
//    
//    RandomStream s;
//    double[] m_ ;
//    
//    m_ = new double[100];
//    
//    xi =  new ChiSquareNoncentralDist (10,20);
//    
//    for(int i=0; i<100; ++i){
//       m_[i] = xi.density(i+1);
//       out.println(m_[i]);
//   }
//                RandomVariateGen gen0 = new RandomVariateGen (new MRG31k3p(), new NormalDist());
//		RandomVariateGen gen1 = new NormalGen (new MRG31k3p());
//		RandomVariateGen gen2 = new NormalGen (new MRG31k3p(), 5.0, 121.4);
//		RandomVariateGen gen3 = new GammaGen (new MRG31k3p(), 2.0, 10.0);
//		RandomVariateGenInt gen4 = new PoissonGen (new MRG31k3p(), 10.0);
//		RandomVariateGen gen5 = new ChiSquareNoncentralGamGen (new LFSR113(), 10.0, 20.0);
//		System.out.println ("Some normal, gamma, , Poisson , ChisquaredNoncenter variates \n");
//		generate(gen0, 3);  // Generate and prints 3 standard normal variates from gen0
//		generate(gen1, 5);  // then 5 more standard normals from gen1
//		generate(gen2, 3);  // then 3 normal variates with mean 5 and standard deviation 121.4 
//		generate(gen3, 2);  // then 2 gamma variates with parameters (2, 10).
//                generate(gen5, 2);  // then 2 gamma variates with parameters (2, 10).
//		int[] arrayP = gen4.nextArrayOfInt (10000); // Generate an array of 10000 Poisson variates of mean 10.
//		Arrays.sort(arrayP);  
//		System.out.println ("50% quantile from Poisson(10) variates: " + arrayP[5000]);		
//		System.out.println ("90% quantile from Poisson(10) variates: " + arrayP[9000]);		
//		System.out.println ("99% quantile from Poisson(10) variates: " + arrayP[9900]);		
//		System.out.println ("99.9% quantile from Poisson(10) variates: " + arrayP[9990]);
//    
//    }
}
    

