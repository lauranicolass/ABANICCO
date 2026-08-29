package abanicco;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JPanel;

/**
 * Interactive polar description panel:
 *   - outer fan ring = categories (click -> % of category in the image)
 *   - inner bubbles  = detected shades, sized by % of image
 *     (click -> radius, angle, L, % and mixing name)
 *   - angle labels outside the ring at every category boundary
 *   - radius scale bar below the wheel
 *   - bubble/scatter display modes (scatter = fixed-size dots, for dense clusters)
 */
public class PolarPanel extends JPanel {

    public interface InfoListener { void onInfo(String[] lines); }

    private final AbaniccoEngine engine;
    private InfoListener listener;
    private boolean scatterMode = false;
    private final List<Shape> wedgeShapes = new ArrayList<>();
    private final List<AbaniccoEngine.Region> wedgeRegions = new ArrayList<>();
    private final List<Ellipse2D> bubbleShapes = new ArrayList<>();
    private final List<AbaniccoEngine.Shade> bubbleShades = new ArrayList<>();
    private static final double RMAX = 150;
    private static final int LABEL_MARGIN = 26;   // room for angle labels
    private static final int SCALEBAR_H = 34;     // room for the radius bar

    public PolarPanel(AbaniccoEngine engine) {
        this.engine = engine;
        setBackground(Color.WHITE);
        setPreferredSize(new Dimension(460, 460 + SCALEBAR_H));
        addMouseListener(new MouseAdapter() {
            @Override public void mouseClicked(MouseEvent e) { hit(e.getPoint()); }
        });
        setToolTipText("");
    }

    /** Scatter mode = all shades drawn as small fixed-size dots. */
    public void setScatterMode(boolean on) {
        if (scatterMode != on) { scatterMode = on; repaint(); }
    }

    public boolean isScatterMode() { return scatterMode; }

    public void setInfoListener(InfoListener l) { listener = l; }

    private void hit(Point p) {
        if (listener == null) return;
        for (int i = bubbleShapes.size() - 1; i >= 0; i--) {
            if (bubbleShapes.get(i).contains(p)) {
                AbaniccoEngine.Shade s = bubbleShades.get(i);
                listener.onInfo(new String[]{
                        "Category: " + s.region,
                        String.format("radius = %.2f | angle = %.2f\u00ba | L = %.2f",
                                s.rho, s.deg, s.L),
                        String.format("%.2f%% of the image", s.percent),
                        s.shadeName});
                return;
            }
        }
        for (int i = 0; i < wedgeShapes.size(); i++) {
            if (wedgeShapes.get(i).contains(p)) {
                AbaniccoEngine.Region r = wedgeRegions.get(i);
                listener.onInfo(new String[]{String.format(
                        "%s: %.2f%% of the image", r.name, r.maskPercentage)});
                return;
            }
        }
    }

    @Override public String getToolTipText(MouseEvent e) {
        Point p = e.getPoint();
        for (int i = bubbleShapes.size() - 1; i >= 0; i--)
            if (bubbleShapes.get(i).contains(p))
                return bubbleShades.get(i).shadeName;
        for (int i = 0; i < wedgeShapes.size(); i++)
            if (wedgeShapes.get(i).contains(p))
                return wedgeRegions.get(i).name;
        return null;
    }

    /** polar (theta rad, rho) -> screen. Angle convention as MATLAB polaraxes. */
    private Point2D toXY(double theta, double rho, double cx, double cy, double scale) {
        double x = cx + Math.cos(theta) * rho * scale;
        double y = cy - Math.sin(theta) * rho * scale;
        return new Point2D.Double(x, y);
    }

    @Override protected void paintComponent(Graphics g0) {
        super.paintComponent(g0);
        Graphics2D g = (Graphics2D) g0;
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        wedgeShapes.clear(); wedgeRegions.clear();
        bubbleShapes.clear(); bubbleShades.clear();

        int size = Math.min(getWidth(), getHeight() - SCALEBAR_H);
        double cx = getWidth() / 2.0;
        double cy = (getHeight() - SCALEBAR_H) / 2.0;
        double scale = (size / 2.0 - LABEL_MARGIN) / RMAX;

        // grid circles
        g.setColor(new Color(230, 230, 230));
        for (int rr = 30; rr <= 150; rr += 30) {
            double rad = rr * scale;
            g.draw(new Ellipse2D.Double(cx - rad, cy - rad, 2 * rad, 2 * rad));
        }

        // category fan ring
        for (AbaniccoEngine.Region reg : engine.regions) {
            double t1 = Math.min(reg.theta1, reg.theta2);
            double t2 = Math.max(reg.theta1, reg.theta2);
            Path2D path = new Path2D.Double();
            int steps = 30;
            for (int i = 0; i <= steps; i++) {
                Point2D p = toXY(t1 + (t2 - t1) * i / steps, reg.rho1, cx, cy, scale);
                if (i == 0) path.moveTo(p.getX(), p.getY());
                else path.lineTo(p.getX(), p.getY());
            }
            for (int i = steps; i >= 0; i--) {
                Point2D p = toXY(t1 + (t2 - t1) * i / steps, reg.rho2, cx, cy, scale);
                path.lineTo(p.getX(), p.getY());
            }
            path.closePath();
            g.setColor(Color.decode(reg.code));
            g.fill(path);
            wedgeShapes.add(path);
            wedgeRegions.add(reg);
        }

        // angle labels at the boundaries between fan categories, outside the ring
        g.setFont(getFont().deriveFont(10f));
        FontMetrics fm = g.getFontMetrics();
        List<Double> boundaries = new ArrayList<>();
        for (AbaniccoEngine.Region reg : engine.regions) {
            if (reg.rho2 < 140) continue;              // only the fan ring
            for (double t : new double[]{reg.theta1, reg.theta2}) {
                double deg = Math.toDegrees(t) % 360;
                if (deg < 0) deg += 360;
                boolean dup = false;
                for (double b : boundaries)
                    if (Math.abs(b - deg) < 0.5 || Math.abs(b - deg) > 359.5) { dup = true; break; }
                if (!dup) boundaries.add(deg);
            }
        }
        g.setColor(new Color(60, 60, 60));
        for (double deg : boundaries) {
            double t = Math.toRadians(deg);
            Point2D p = toXY(t, RMAX + 9, cx, cy, scale);
            String txt = String.format("%.0f\u00ba", deg);
            double tw = fm.stringWidth(txt);
            // shift the anchor so text stays outside the wheel on every side
            double x = p.getX() - tw / 2 + Math.cos(t) * tw / 2;
            double y = p.getY() + fm.getAscent() / 2.0 - Math.sin(t) * fm.getAscent() / 2.0;
            g.drawString(txt, (float) x, (float) y);
        }

        // guide lines
        g.setColor(Color.BLACK);
        g.setStroke(new BasicStroke(0.7f));
        for (double[] ln : engine.lines) {
            Point2D a = toXY(ln[0], ln[2] * 98, cx, cy, scale);
            Point2D b = toXY(ln[1], ln[3] * 98, cx, cy, scale);
            g.draw(new Line2D.Double(a, b));
        }
        double r3 = 3 * scale;
        g.draw(new Ellipse2D.Double(cx - r3, cy - r3, 2 * r3, 2 * r3));
        // brown arc at rho 35
        AbaniccoEngine.Region brown = engine.regions.get(engine.regions.size() - 2);
        Path2D arc = new Path2D.Double();
        for (int i = 0; i <= 40; i++) {
            double t = brown.theta1 + (brown.theta2 - brown.theta1) * i / 40.0;
            Point2D p = toXY(t, 35, cx, cy, scale);
            if (i == 0) arc.moveTo(p.getX(), p.getY());
            else arc.lineTo(p.getX(), p.getY());
        }
        g.draw(arc);
        // base stars (small black diamonds)
        for (AbaniccoEngine.Region reg : engine.regions) {
            Point2D p = toXY(reg.base, reg.rho1, cx, cy, scale);
            g.fill(new Ellipse2D.Double(p.getX() - 2.5, p.getY() - 2.5, 5, 5));
        }

        // shades
        double maxPer = 1e-9;
        for (AbaniccoEngine.Shade s : engine.allShades())
            maxPer = Math.max(maxPer, s.percent);
        for (AbaniccoEngine.Shade s : engine.allShades()) {
            double d = scatterMode ? 7 : 8 + 26 * s.percent / maxPer;
            Point2D p = toXY(s.theta, s.rho, cx, cy, scale);
            Ellipse2D e = new Ellipse2D.Double(p.getX() - d / 2, p.getY() - d / 2, d, d);
            g.setColor(new Color(s.r, s.g, s.b));
            g.fill(e);
            g.setColor(Color.BLACK);
            g.setStroke(new BasicStroke(1f));
            g.draw(e);
            bubbleShapes.add(e);
            bubbleShades.add(s);
        }

        // radius scale bar below the wheel (same pixel scale as the plot)
        double barY = getHeight() - SCALEBAR_H + 10;
        double barX0 = cx - RMAX * scale / 2.0;   // centered, spans rho 0..150
        g.setColor(Color.BLACK);
        g.setStroke(new BasicStroke(1f));
        g.draw(new Line2D.Double(barX0, barY, barX0 + RMAX * scale, barY));
        g.setFont(getFont().deriveFont(10f));
        FontMetrics fm2 = g.getFontMetrics();
        for (int rr = 0; rr <= 150; rr += 30) {
            double x = barX0 + rr * scale;
            g.draw(new Line2D.Double(x, barY - 3, x, barY + 3));
            String txt = Integer.toString(rr);
            g.drawString(txt, (float) (x - fm2.stringWidth(txt) / 2.0),
                    (float) (barY + 5 + fm2.getAscent()));
        }
        String cap = "radius \u03c1 = \u221a(a\u00b2+b\u00b2)";
        g.drawString(cap, (float) (barX0 + RMAX * scale + 10), (float) (barY + 4));
    }
}
