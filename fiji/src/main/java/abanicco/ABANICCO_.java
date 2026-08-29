package abanicco;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * ABANICCO for FIJI / ImageJ.
 *
 * Workflow (per the "Road to FIJI" design):
 *   1. Input: current image + dialog for #colors per category and morphology
 *   2. General visualization: original image, category segmentations,
 *      interactive polar description
 *   3. Post-processing: change boundaries, create clusters,
 *      change #colors, apply morphology
 *   4. Export checklist
 */
public class ABANICCO_ implements PlugIn {

    private AbaniccoEngine engine;
    private ImagePlus imp;
    private int[] pixels;
    private int w, h;
    private int numShades = 15, morpho = 10;
    private final List<String> log = new ArrayList<>();

    private JFrame frame;
    private JPanel segGrid;
    private JTextArea infoArea;
    private PolarPanel polar;

    @Override
    public void run(String arg) {
        imp = IJ.getImage();
        if (imp == null) {
            IJ.error("ABANICCO", "Open an image first.");
            return;
        }
        // ---- 1. Input dialog ------------------------------------------------
        GenericDialog gd = new GenericDialog("ABANICCO - Input");
        gd.addNumericField("Colors per category (bubbles):", numShades, 0);
        gd.addNumericField("Morphology (px, 10 for ~500x500 images):", morpho, 0);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        numShades = (int) gd.getNextNumber();
        morpho = (int) gd.getNextNumber();

        try {
            engine = AbaniccoEngine.loadDefault();
        } catch (Exception e) {
            IJ.handleException(e);
            return;
        }

        ColorProcessor cp = (ColorProcessor) imp.getProcessor().convertToRGB();
        w = cp.getWidth();
        h = cp.getHeight();
        pixels = (int[]) cp.getPixels();

        IJ.showStatus("ABANICCO: analyzing...");
        engine.analyze(pixels, w, h, numShades, morpho);
        log.add("analyze(numShades=" + numShades + ", morpho=" + morpho + ")");
        IJ.showStatus("ABANICCO: done");

        SwingUtilities.invokeLater(this::buildUI);
    }

    // ------------------------------------------------------------------- UI
    private void buildUI() {
        frame = new JFrame("ABANICCO - " + imp.getTitle());
        frame.setLayout(new BorderLayout(6, 6));
        frame.getContentPane().setBackground(Color.WHITE);

        // LEFT: original image
        JPanel left = new JPanel(new BorderLayout());
        left.setBorder(new TitledBorder("Original image"));
        left.add(new JLabel(new ImageIcon(scaled(toBuffered(pixels, w, h), 260))),
                BorderLayout.CENTER);
        frame.add(left, BorderLayout.WEST);

        // CENTER: segmentation grid (click a tile -> distribution + shade maps)
        segGrid = new JPanel(new GridLayout(0, 4, 4, 4));
        segGrid.setBackground(Color.WHITE);
        refreshSegGrid();
        JScrollPane sp = new JScrollPane(segGrid);
        sp.setBorder(new TitledBorder("Category segmentations (click to inspect)"));
        frame.add(sp, BorderLayout.CENTER);

        // RIGHT: polar description + info
        JPanel right = new JPanel(new BorderLayout());
        polar = new PolarPanel(engine);
        polar.setBorder(new TitledBorder("Polar description (click fan or bubbles)"));
        infoArea = new JTextArea(5, 30);
        infoArea.setEditable(false);
        infoArea.setText("Click a shade or a category in the fan.");
        polar.setInfoListener(lines -> infoArea.setText(String.join("\n", lines)));
        JPanel modeBar = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 2));
        modeBar.setBackground(Color.WHITE);
        modeBar.add(new JLabel("Shades:"));
        JRadioButton rbBubble = new JRadioButton("Bubble plot (size = %)", true);
        JRadioButton rbScatter = new JRadioButton("Scatter plot (fixed size)");
        rbBubble.setBackground(Color.WHITE);
        rbScatter.setBackground(Color.WHITE);
        ButtonGroup grp = new ButtonGroup();
        grp.add(rbBubble); grp.add(rbScatter);
        rbBubble.addActionListener(e -> polar.setScatterMode(false));
        rbScatter.addActionListener(e -> polar.setScatterMode(true));
        modeBar.add(rbBubble); modeBar.add(rbScatter);
        right.add(modeBar, BorderLayout.NORTH);
        right.add(polar, BorderLayout.CENTER);
        right.add(new JScrollPane(infoArea), BorderLayout.SOUTH);
        frame.add(right, BorderLayout.EAST);

        // BOTTOM: post-processing + export
        JPanel bottom = new JPanel(new FlowLayout(FlowLayout.LEFT));
        bottom.setBackground(Color.WHITE);
        JButton bBounds = new JButton("Change boundaries...");
        JButton bCluster = new JButton("Create clusters...");
        JButton bColors = new JButton("Change #colors / morphology...");
        JButton bExport = new JButton("Export...");
        bBounds.addActionListener(e -> changeBoundaries());
        bCluster.addActionListener(e -> createClusters());
        bColors.addActionListener(e -> rerunDialog());
        bExport.addActionListener(e -> exportDialog());
        bottom.add(bBounds);
        bottom.add(bCluster);
        bottom.add(bColors);
        bottom.add(bExport);
        frame.add(bottom, BorderLayout.SOUTH);

        frame.pack();
        frame.setSize(Math.min(1380, frame.getWidth()), 760);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private void refreshSegGrid() {
        segGrid.removeAll();
        for (AbaniccoEngine.Region reg : engine.regions) {
            BufferedImage seg = toBuffered(maskPixels(reg.mask), w, h);
            JLabel lab = new JLabel(String.format("<html><center>%s<br>(%.1f%%)</center></html>",
                    reg.name, reg.maskPercentage),
                    new ImageIcon(scaled(seg, 150)), SwingConstants.CENTER);
            lab.setVerticalTextPosition(SwingConstants.BOTTOM);
            lab.setHorizontalTextPosition(SwingConstants.CENTER);
            lab.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            lab.addMouseListener(new MouseAdapter() {
                @Override public void mouseClicked(MouseEvent e) { inspectCategory(reg); }
            });
            segGrid.add(lab);
        }
        segGrid.revalidate();
        segGrid.repaint();
    }

    /** Click on a segmentation: show it in ImageJ + individual color maps. */
    private void inspectCategory(AbaniccoEngine.Region reg) {
        infoArea.setText(String.format("%s: %.2f%% of the image",
                reg.name, reg.maskPercentage));
        new ImagePlus("Segmented " + reg.name,
                new ColorProcessor(toBuffered(maskPixels(reg.mask), w, h))).show();
        for (int i = 0; i < reg.shades.size(); i++) {
            AbaniccoEngine.Shade s = reg.shades.get(i);
            ByteProcessor bp = new ByteProcessor(w, h);
            for (int p = 0; p < s.mask.length; p++)
                if (s.mask[p]) bp.set(p % w, p / w, 255);
            new ImagePlus(String.format("%s shade %d - %s",
                    reg.name, i + 1, s.shadeName), bp).show();
        }
    }

    // ------------------------------------------------------ post-processing
    private void changeBoundaries() {
        GenericDialog gd = new GenericDialog("ABANICCO - Change boundaries");
        gd.addSlider("Theta 1 (deg)", 0, 360, 90);
        gd.addSlider("Theta 2 (deg)", 0, 360, 180);
        gd.addSlider("Rho min", 0, 130, 30);
        gd.addSlider("Rho max", 0, 200, 200);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        double t1 = gd.getNextNumber(), t2 = gd.getNextNumber();
        double r1 = gd.getNextNumber(), r2 = gd.getNextNumber();
        boolean[] m = AbaniccoEngine.customRegion(pixels, t1, t2, r1, r2);
        new ImagePlus(String.format("Custom region theta(%.0f-%.0f) rho(%.0f-%.0f)",
                t1, t2, r1, r2),
                new ColorProcessor(toBuffered(maskPixels(m), w, h))).show();
        log.add(String.format("custom region: theta(%.0f,%.0f) rho(%.0f,%.0f)",
                t1, t2, r1, r2));
    }

    private void createClusters() {
        GenericDialog gd = new GenericDialog("ABANICCO - Create cluster");
        for (AbaniccoEngine.Region r : engine.regions)
            gd.addCheckbox(r.name, false);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        List<String> sel = new ArrayList<>();
        for (AbaniccoEngine.Region r : engine.regions)
            if (gd.getNextBoolean()) sel.add(r.name);
        if (sel.isEmpty()) return;
        boolean[] m = engine.cluster(sel.toArray(new String[0]), w * h);
        new ImagePlus("Cluster: " + String.join("+", sel),
                new ColorProcessor(toBuffered(maskPixels(m), w, h))).show();
        log.add("cluster: " + String.join("+", sel));
    }

    private void rerunDialog() {
        GenericDialog gd = new GenericDialog("ABANICCO - Re-run");
        gd.addNumericField("Colors per category:", numShades, 0);
        gd.addNumericField("Morphology (px):", morpho, 0);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        numShades = (int) gd.getNextNumber();
        morpho = (int) gd.getNextNumber();
        IJ.showStatus("ABANICCO: re-analyzing...");
        engine.analyze(pixels, w, h, numShades, morpho);
        log.add("analyze(numShades=" + numShades + ", morpho=" + morpho + ")");
        refreshSegGrid();
        polar.repaint();
        IJ.showStatus("ABANICCO: done");
    }

    // --------------------------------------------------------------- export
    private void exportDialog() {
        GenericDialog gd = new GenericDialog("ABANICCO - Export");
        gd.addCheckbox("Polar description (PNG)", true);
        gd.addCheckbox("Segmentations of each category (PNG)", true);
        gd.addCheckbox("Individual color masks of each category (PNG)", false);
        gd.addCheckbox("Table of all detected colors (CSV)", true);
        gd.addCheckbox("List of modifications (TXT)", true);
        gd.addCheckbox("Code for MATLAB / Python (batch)", true);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        boolean ePolar = gd.getNextBoolean(), eSeg = gd.getNextBoolean(),
                eInd = gd.getNextBoolean(), eTab = gd.getNextBoolean(),
                eLog = gd.getNextBoolean(), eCode = gd.getNextBoolean();
        String dir = IJ.getDirectory("Choose export folder");
        if (dir == null) return;
        try {
            if (ePolar) {
                BufferedImage bi = new BufferedImage(800, 800,
                        BufferedImage.TYPE_INT_RGB);
                Graphics2D g = bi.createGraphics();
                g.setColor(Color.WHITE);
                g.fillRect(0, 0, 800, 800);
                polar.setSize(800, 800);
                polar.paint(g);
                g.dispose();
                javax.imageio.ImageIO.write(bi, "png",
                        new File(dir, "polar_description.png"));
            }
            if (eSeg)
                for (AbaniccoEngine.Region r : engine.regions)
                    javax.imageio.ImageIO.write(
                            toBuffered(maskPixels(r.mask), w, h), "png",
                            new File(dir, "segmentation_"
                                    + r.name.replace(' ', '_') + ".png"));
            if (eInd)
                for (AbaniccoEngine.Region r : engine.regions)
                    for (int i = 0; i < r.shades.size(); i++) {
                        BufferedImage bm = new BufferedImage(w, h,
                                BufferedImage.TYPE_BYTE_GRAY);
                        boolean[] m = r.shades.get(i).mask;
                        for (int p = 0; p < m.length; p++)
                            if (m[p]) bm.setRGB(p % w, p / w, 0xffffff);
                        javax.imageio.ImageIO.write(bm, "png", new File(dir,
                                r.name.replace(' ', '_') + "_shade" + (i + 1) + ".png"));
                    }
            if (eTab) {
                ResultsTable rt = new ResultsTable();
                for (AbaniccoEngine.Shade s : engine.allShades()) {
                    rt.incrementCounter();
                    rt.addValue("Category", s.region);
                    rt.addValue("ShadeName", s.shadeName);
                    rt.addValue("Radius", s.rho);
                    rt.addValue("Angle_deg", s.deg);
                    rt.addValue("L", s.L);
                    rt.addValue("Percent_of_image", s.percent);
                    rt.addValue("Hex", String.format("#%02x%02x%02x",
                            Math.round(s.r * 255), Math.round(s.g * 255),
                            Math.round(s.b * 255)));
                }
                rt.save(new File(dir, "shades.csv").getAbsolutePath());
                rt.show("ABANICCO shades");
            }
            if (eLog)
                try (PrintWriter pw = new PrintWriter(
                        new FileWriter(new File(dir, "modifications.txt")))) {
                    for (String s : log) pw.println(s);
                }
            if (eCode) {
                try (PrintWriter pw = new PrintWriter(
                        new FileWriter(new File(dir, "reproduce.py")))) {
                    pw.println("# Auto-generated by ABANICCO for FIJI - batch ready");
                    pw.println("import numpy as np");
                    pw.println("from PIL import Image");
                    pw.println("from abanicco import analyze, shades_to_records");
                    pw.println("import pandas as pd");
                    pw.println();
                    pw.println("im = np.array(Image.open(r'" + imp.getTitle()
                            + "').convert('RGB'))");
                    pw.println("regions, rows = analyze(im, numshades="
                            + numShades + ", morpho=" + morpho + ")");
                    pw.println("pd.DataFrame(shades_to_records(rows))"
                            + ".to_csv('shades.csv', index=False)");
                }
                try (PrintWriter pw = new PrintWriter(
                        new FileWriter(new File(dir, "reproduce.m")))) {
                    pw.println("% Auto-generated by ABANICCO for FIJI - batch ready");
                    pw.println("load('Regions_Angles0.mat'); load('Regions_Angles1.mat');");
                    pw.println("an_im = imread('" + imp.getTitle() + "');");
                    pw.println("[Results,~,subTables] = AnalyzeColorImage_final2_short("
                            + "an_im, Regions_Angles, Regions_Angles1, "
                            + numShades + ", " + morpho + ", 1);");
                }
            }
            IJ.showMessage("ABANICCO", "Exported to " + dir);
        } catch (Exception e) {
            IJ.handleException(e);
        }
    }

    // -------------------------------------------------------------- helpers
    private int[] maskPixels(boolean[] mask) {
        int[] out = new int[pixels.length];
        for (int i = 0; i < pixels.length; i++)
            out[i] = mask[i] ? pixels[i] : 0;
        return out;
    }

    private static BufferedImage toBuffered(int[] px, int w, int h) {
        BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        bi.setRGB(0, 0, w, h, px, 0, w);
        return bi;
    }

    private static Image scaled(BufferedImage bi, int maxDim) {
        double f = (double) maxDim / Math.max(bi.getWidth(), bi.getHeight());
        return bi.getScaledInstance((int) (bi.getWidth() * f),
                (int) (bi.getHeight() * f), Image.SCALE_SMOOTH);
    }
}
