package abanicco;

import java.util.ArrayDeque;

/** bwareaopen equivalent: remove 8-connected components smaller than minSize. */
public class Morphology {

    public static boolean[] areaOpen(boolean[] mask, int w, int h, int minSize) {
        if (minSize <= 1) return mask;
        boolean[] out = new boolean[mask.length];
        boolean[] visited = new boolean[mask.length];
        int[] stackBuf = new int[mask.length];
        ArrayDeque<Integer> comp = new ArrayDeque<>();
        for (int start = 0; start < mask.length; start++) {
            if (!mask[start] || visited[start]) continue;
            // flood fill (8-connectivity)
            int top = 0;
            stackBuf[top++] = start;
            visited[start] = true;
            comp.clear();
            while (top > 0) {
                int p = stackBuf[--top];
                comp.add(p);
                int y = p / w, x = p % w;
                for (int dy = -1; dy <= 1; dy++) {
                    int ny = y + dy;
                    if (ny < 0 || ny >= h) continue;
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dx == 0 && dy == 0) continue;
                        int nx = x + dx;
                        if (nx < 0 || nx >= w) continue;
                        int q = ny * w + nx;
                        if (mask[q] && !visited[q]) {
                            visited[q] = true;
                            stackBuf[top++] = q;
                        }
                    }
                }
            }
            if (comp.size() >= minSize)
                for (int p : comp) out[p] = true;
        }
        return out;
    }
}
