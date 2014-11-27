using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoMath
{
    public static class LinearEstimator
    {
        public static double EstimateVariance(Func<double, double> Cov, Func<double, double, double, double, double> Dist, double sill, double targetLat, double targetLon, int[] nodeIndeces, double[] nodeWeights, double[] nodeLats, double[] nodeLons)
        {
            //var = cov(0)+ sum sum (w[i]*w[j]*cov(i,j))-2.0*sum(w[i]*cov(x,i))
            double cov_at_0 = sill;

            double acc = cov_at_0; //cov(0)                    
            for (int i = 0; i < nodeWeights.Length; i++)
            {
                double w = nodeWeights[i];
                int idx1 = nodeIndeces[i];
                double lat1 = nodeLats[idx1];
                double lon1 = nodeLons[idx1];
                for (int j = 0; j < i; j++)
                {
                    int idx2 = nodeIndeces[j];
                    double lat2 = nodeLats[idx2];
                    double lon2 = nodeLons[idx2];
                    double dist = Dist(lat1, lon1, lat2, lon2);
                    double cov = Cov(dist);
                    acc += 2.0 * w * nodeWeights[j] * cov;
                }
                acc += w * w * cov_at_0; //diagonal elements
                double dist2 = Dist(lat1, lon1, targetLat, targetLon);
                double cov2 = Cov(dist2);
                acc -= 2.0 * w * cov2;
            }
            return acc;
        }
    }
}
