package abanicco;

import java.util.List;
import java.util.Random;

/** Minimal k-means (k-means++ init) for RGB quantization (rgb2ind analog). */
public class KMeans {

    /** @return k (or fewer) cluster centers as float[3] RGB in 0..1. */
    public static float[][] fit(List<float[]> pts, int k, int iters, long seed) {
        int n = pts.size();
        if (n == 0) return new float[0][];
        k = Math.min(k, n);
        Random rnd = new Random(seed);
        float[][] c = new float[k][3];

        // k-means++ init
        c[0] = pts.get(rnd.nextInt(n)).clone();
        double[] d2 = new double[n];
        for (int ci = 1; ci < k; ci++) {
            double sum = 0;
            for (int i = 0; i < n; i++) {
                double best = Double.MAX_VALUE;
                float[] p = pts.get(i);
                for (int j = 0; j < ci; j++) best = Math.min(best, dist2(p, c[j]));
                d2[i] = best;
                sum += best;
            }
            double target = rnd.nextDouble() * sum, acc = 0;
            int pick = n - 1;
            for (int i = 0; i < n; i++) {
                acc += d2[i];
                if (acc >= target) { pick = i; break; }
            }
            c[ci] = pts.get(pick).clone();
        }

        int[] assign = new int[n];
        for (int it = 0; it < iters; it++) {
            boolean changed = false;
            for (int i = 0; i < n; i++) {
                float[] p = pts.get(i);
                int best = 0;
                double bd = Double.MAX_VALUE;
                for (int j = 0; j < k; j++) {
                    double d = dist2(p, c[j]);
                    if (d < bd) { bd = d; best = j; }
                }
                if (assign[i] != best) { assign[i] = best; changed = true; }
            }
            float[][] sum = new float[k][3];
            int[] cnt = new int[k];
            for (int i = 0; i < n; i++) {
                float[] p = pts.get(i);
                int a = assign[i];
                sum[a][0] += p[0]; sum[a][1] += p[1]; sum[a][2] += p[2];
                cnt[a]++;
            }
            for (int j = 0; j < k; j++)
                if (cnt[j] > 0)
                    for (int d = 0; d < 3; d++) c[j][d] = sum[j][d] / cnt[j];
            if (!changed) break;
        }
        // drop empty clusters
        int[] cnt = new int[k];
        for (int i = 0; i < n; i++) cnt[assign[i]]++;
        int live = 0;
        for (int j = 0; j < k; j++) if (cnt[j] > 0) live++;
        float[][] out = new float[live][];
        int o = 0;
        for (int j = 0; j < k; j++) if (cnt[j] > 0) out[o++] = c[j];
        return out;
    }

    private static double dist2(float[] a, float[] b) {
        double dr = a[0] - b[0], dg = a[1] - b[1], db = a[2] - b[2];
        return dr * dr + dg * dg + db * db;
    }
}
