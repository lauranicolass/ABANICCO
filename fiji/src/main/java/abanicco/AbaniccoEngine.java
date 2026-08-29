package abanicco;

import java.io.*;
import java.util.*;

/**
 * ABANICCO core engine: pure computation, no UI.
 * Port of AnalyzeColorImage_final2_short.m
 * Reference: https://www.mdpi.com/1424-8220/23/6/3338
 */
public class AbaniccoEngine {

    // ------------------------------------------------------------- model
    public static class Region {
        public String name, code;
        public double theta1, theta2, rho1, rho2, base;
        // results
        public boolean[] mask;            // h*w
        public double maskPercentage;
        public List<Shade> shades = new ArrayList<>();
    }

    public static class Shade {
        public String region, code, shadeName = "";
        public double theta, rho, L, deg, percent; // percent of whole image
        public float r, g, b;                       // 0..1
        public boolean[] mask;
    }

    public static class Interval {
        public String color;
        public double p1, p2, x1, x2;
    }

    public List<Region> regions = new ArrayList<>();
    public List<Interval> intervals = new ArrayList<>();
    public List<Interval> intervalsRad = new ArrayList<>();
    public double[][] lines;   // [n][4]: t1,t2,r1,r2 guide lines
    private short[] lut;       // 256*256, bit i = membership of region i

    // ------------------------------------------------------ data loading
    public static AbaniccoEngine loadDefault() throws IOException {
        AbaniccoEngine e = new AbaniccoEngine();
        e.loadResources();
        return e;
    }

    private BufferedReader res(String name) {
        InputStream in = getClass().getResourceAsStream("/abanicco/" + name);
        if (in == null)
            throw new RuntimeException("Missing resource /abanicco/" + name);
        return new BufferedReader(new InputStreamReader(in));
    }

    private void loadResources() throws IOException {
        try (BufferedReader br = res("regions.csv")) {
            br.readLine();
            String s;
            while ((s = br.readLine()) != null) {
                String[] t = s.split(",");
                Region r = new Region();
                r.name = t[0]; r.code = t[1];
                r.theta1 = Double.parseDouble(t[2]);
                r.theta2 = Double.parseDouble(t[3]);
                r.rho1 = Double.parseDouble(t[4]);
                r.rho2 = Double.parseDouble(t[5]);
                r.base = Double.parseDouble(t[6]);
                regions.add(r);
            }
        }
        lut = new short[256 * 256];
        try (DataInputStream in = new DataInputStream(new BufferedInputStream(
                getClass().getResourceAsStream("/abanicco/region_lut.bin")))) {
            for (int i = 0; i < lut.length; i++) lut[i] = in.readShort();
        }
        intervals = readIntervals("intervals.csv");
        intervalsRad = readIntervals("intervals_rad.csv");
        List<double[]> ls = new ArrayList<>();
        try (BufferedReader br = res("lines.csv")) {
            br.readLine();
            String s;
            while ((s = br.readLine()) != null) {
                String[] t = s.split(",");
                ls.add(new double[]{Double.parseDouble(t[0]), Double.parseDouble(t[1]),
                        Double.parseDouble(t[2]), Double.parseDouble(t[3])});
            }
        }
        lines = ls.toArray(new double[0][]);
    }

    private List<Interval> readIntervals(String file) throws IOException {
        List<Interval> out = new ArrayList<>();
        try (BufferedReader br = res(file)) {
            br.readLine();
            String s;
            while ((s = br.readLine()) != null) {
                String[] t = s.split(",");
                Interval iv = new Interval();
                iv.color = t[0];
                iv.p1 = Double.parseDouble(t[1]);
                iv.p2 = Double.parseDouble(t[2]);
                iv.x1 = Double.parseDouble(t[3]);
                iv.x2 = Double.parseDouble(t[4]);
                out.add(iv);
            }
        }
        return out;
    }

    // ---------------------------------------------------- color conversion
    /** sRGB (0..255) -> CIELAB (D65), matching MATLAB rgb2lab. */
    public static void rgb2lab(int rgb, double[] lab) {
        double r = srgbInv(((rgb >> 16) & 0xff) / 255.0);
        double g = srgbInv(((rgb >> 8) & 0xff) / 255.0);
        double b = srgbInv((rgb & 0xff) / 255.0);
        double X = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b;
        double Y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
        double Z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b;
        double fx = labF(X / 0.95047), fy = labF(Y), fz = labF(Z / 1.08883);
        lab[0] = 116 * fy - 16;
        lab[1] = 500 * (fx - fy);
        lab[2] = 200 * (fy - fz);
    }

    private static double srgbInv(double c) {
        return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
    }

    private static double labF(double t) {
        return t > 0.008856 ? Math.cbrt(t) : (7.787 * t + 16.0 / 116.0);
    }

    public static double wrapTo360(double deg) {
        double d = deg % 360.0;
        return d < 0 ? d + 360.0 : d;
    }

    // ------------------------------------------------------------ analysis
    /**
     * Full analysis. pixels: packed RGB ints (0xRRGGBB), row-major.
     * numShades: colors detected per category; morpho: min object size.
     */
    public void analyze(int[] pixels, int w, int h, int numShades, int morpho) {
        int n = w * h;
        // per-pixel LAB a,b -> LUT index
        short[] member = new short[n];
        double[] lab = new double[3];
        for (int i = 0; i < n; i++) {
            rgb2lab(pixels[i], lab);
            int ai = (int) Math.floor(lab[1] + 127);
            int bi = (int) Math.floor(lab[2] + 127);
            if (ai < 0) ai = 0; if (ai > 255) ai = 255;
            if (bi < 0) bi = 0; if (bi > 255) bi = 255;
            member[i] = lut[ai * 256 + bi];
        }
        for (int ri = 0; ri < regions.size(); ri++) {
            Region reg = regions.get(ri);
            boolean[] m = new boolean[n];
            int bit = 1 << ri;
            for (int i = 0; i < n; i++) m[i] = (member[i] & bit) != 0;
            m = Morphology.areaOpen(m, w, h, morpho);
            boolean[] inv = new boolean[n];
            for (int i = 0; i < n; i++) inv[i] = !m[i];
            inv = Morphology.areaOpen(inv, w, h, 2 * morpho);
            for (int i = 0; i < n; i++) m[i] = !inv[i];
            reg.mask = m;
            int cnt = 0;
            for (boolean v : m) if (v) cnt++;
            reg.maskPercentage = 100.0 * cnt / n;
            reg.shades = computeShades(reg, pixels, n, numShades, cnt);
        }
        nameShades();
    }

    private List<Shade> computeShades(Region reg, int[] pixels, int n,
                                      int numShades, int regionPixels) {
        List<Shade> out = new ArrayList<>();
        if (regionPixels == 0) return out;
        // gather region pixels (subsample to cap k-means cost)
        int step = Math.max(1, regionPixels / 60000);
        List<float[]> sample = new ArrayList<>();
        int seen = 0;
        for (int i = 0; i < n; i++) {
            if (!reg.mask[i]) continue;
            if ((seen++ % step) != 0) continue;
            int p = pixels[i];
            sample.add(new float[]{((p >> 16) & 0xff) / 255f,
                    ((p >> 8) & 0xff) / 255f, (p & 0xff) / 255f});
        }
        float[][] centers = KMeans.fit(sample, numShades, 25, 0L);
        // assign all region pixels to the nearest center
        int k = centers.length;
        int[] counts = new int[k];
        int[] assign = new int[n];
        Arrays.fill(assign, -1);
        for (int i = 0; i < n; i++) {
            if (!reg.mask[i]) continue;
            int p = pixels[i];
            float r = ((p >> 16) & 0xff) / 255f, g = ((p >> 8) & 0xff) / 255f,
                    b = (p & 0xff) / 255f;
            int best = 0;
            float bd = Float.MAX_VALUE;
            for (int c = 0; c < k; c++) {
                float dr = r - centers[c][0], dg = g - centers[c][1],
                        db = b - centers[c][2];
                float d = dr * dr + dg * dg + db * db;
                if (d < bd) { bd = d; best = c; }
            }
            assign[i] = best;
            counts[best]++;
        }
        for (int c = 0; c < k; c++) {
            double pctOfRegion = 100.0 * counts[c] / regionPixels;
            if (pctOfRegion <= 5.0) continue;
            Shade s = new Shade();
            s.region = reg.name;
            s.code = reg.code;
            s.r = centers[c][0]; s.g = centers[c][1]; s.b = centers[c][2];
            int rgb = ((int) Math.round(s.r * 255) << 16)
                    | ((int) Math.round(s.g * 255) << 8)
                    | (int) Math.round(s.b * 255);
            double[] lab = new double[3];
            rgb2lab(rgb, lab);
            s.L = lab[0];
            s.theta = Math.atan2(lab[1], lab[2]); // cart2pol(b,a)
            s.rho = Math.hypot(lab[1], lab[2]);
            s.deg = wrapTo360(Math.toDegrees(s.theta));
            s.percent = 100.0 * counts[c] / n;
            s.mask = new boolean[n];
            for (int i = 0; i < n; i++) s.mask[i] = assign[i] == c;
            out.add(s);
        }
        return out;
    }

    // -------------------------------------------------------- shade naming
    private void nameShades() {
        // radial interval pair with duplicated bounds (Brown / Pure)
        int z1 = -1, z2 = -1;
        for (int i = 0; i < intervalsRad.size() && z2 < 0; i++)
            for (int j = i + 1; j < intervalsRad.size(); j++)
                if (intervalsRad.get(i).x1 == intervalsRad.get(j).x1
                        && intervalsRad.get(i).x2 == intervalsRad.get(j).x2) {
                    z1 = i; z2 = j; break;
                }
        for (Region reg : regions)
            for (Shade s : reg.shades)
                s.shadeName = nameOne(s, z1, z2);
    }

    private String nameOne(Shade s, int z1, int z2) {
        double th = s.deg, R = s.rho;
        List<Integer> idxs = new ArrayList<>();
        for (int k = 0; k < intervals.size(); k++)
            if (th < intervals.get(k).x2 && th > intervals.get(k).x1)
                idxs.add(k);
        if (idxs.isEmpty()) return s.region;
        idxs.sort(Collections.reverseOrder());
        if (!intervals.get(idxs.get(0)).color.equals(s.region))
            Collections.sort(idxs);
        Interval iv1 = intervals.get(idxs.get(0));
        double val1 = th * iv1.p1 + iv1.p2;
        Interval iv2 = null;
        double val2 = 0;
        if (idxs.size() > 1) {
            iv2 = intervals.get(idxs.get(1));
            val2 = th * iv2.p1 + iv2.p2;
        }
        boolean withBrown = false;
        double valPure = 0, valBrown = 0;
        if (th > 0 && th < 90 && z1 >= 0) {
            Interval a = intervalsRad.get(z1), b = intervalsRad.get(z2);
            if (R < a.x2 && R > a.x1) {
                withBrown = true;
                double v1 = R * a.p1 + a.p2, v2 = R * b.p1 + b.p2;
                if (a.color.equals("Pure")) { valPure = v1; valBrown = v2; }
                else { valPure = v2; valBrown = v1; }
            }
        }
        if (!withBrown) {
            if (iv2 == null)
                return String.format("%.1f%% of %s", val1 * 100, iv1.color);
            return String.format("%.1f%% of %s and %.1f%% of %s",
                    val1 * 100, iv1.color, val2 * 100, iv2.color);
        }
        if (iv2 == null)
            return String.format("%.1f%% of Brown, and %.1f%% of %s",
                    valBrown * 100, valPure * val1 * 100, iv1.color);
        return String.format("%.1f%% of Brown, %.1f%% of %s, and %.1f%% of %s",
                valBrown * 100, valPure * val1 * 100, iv1.color,
                valPure * val2 * 100, iv2.color);
    }

    // -------------------------------------------------- post-processing
    /** Merge category masks into a cluster mask. */
    public boolean[] cluster(String[] names, int n) {
        boolean[] out = new boolean[n];
        for (String nm : names)
            for (Region r : regions)
                if (r.name.equals(nm) && r.mask != null)
                    for (int i = 0; i < n; i++) out[i] |= r.mask[i];
        return out;
    }

    /** Custom AB-plane region: theta in degrees (wrapped), rho bounds. */
    public static boolean[] customRegion(int[] pixels, double t1, double t2,
                                         double rhoMin, double rhoMax) {
        boolean[] m = new boolean[pixels.length];
        double[] lab = new double[3];
        for (int i = 0; i < pixels.length; i++) {
            rgb2lab(pixels[i], lab);
            double th = wrapTo360(Math.toDegrees(Math.atan2(lab[1], lab[2])));
            double rho = Math.hypot(lab[1], lab[2]);
            boolean inT = t1 <= t2 ? (th > t1 && th < t2) : (th > t1 || th < t2);
            m[i] = inT && rho > rhoMin && rho < rhoMax;
        }
        return m;
    }

    public List<Shade> allShades() {
        List<Shade> out = new ArrayList<>();
        for (Region r : regions) out.addAll(r.shades);
        return out;
    }
}
