import java.io.*;
import java.util.*;
import org.glycoinfo.WURCSFramework.util.map.MAPFactory;
import org.glycoinfo.WURCSFramework.wurcs.map.*;

/** Extracts objective structural features from a WURCS MAP code using the
 *  wurcsframework parser. Input: MAP code <TAB> ... (rest passed through).
 *  Output: TSV of features. */
public class MapFeatures {

    static List<MAPAtomAbstract> V = new ArrayList<>();     // vertices (non-cyclic atoms)
    static Map<MAPAtomAbstract,Integer> IDX = new IdentityHashMap<>();
    static List<int[]> E = new ArrayList<>();               // edges [u,v,bondOrder]
    static List<int[]> CLOSE = new ArrayList<>();           // ring-closure edges

    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));
        out.println(String.join("\t", "map","rest","status","nAtom","nC","nO","nN","nS","nP","nHal",
            "nAromAtom","nRing","nAromRing","ringSizes","nCarbonyl","nCarboxylLike",
            "maxChainC","attach","nStar"));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            int t = line.indexOf('\t');
            String map  = (t < 0) ? line : line.substring(0, t);
            String rest = (t < 0) ? ""   : line.substring(t + 1);
            try {
                out.println(features(map, rest));
            } catch (Throwable e) {
                out.println(map + "\t" + rest + "\tERROR" + "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t");
            }
        }
        out.flush();
    }

    static String features(String map, String rest) throws Exception {
        V.clear(); IDX.clear(); E.clear(); CLOSE.clear();
        MAPGraph g = new MAPFactory(map).getMAPGraph();
        collect(g);

        // --- build adjacency ---
        int n = V.size();
        List<List<Integer>> adj = new ArrayList<>();
        for (int i = 0; i < n; i++) adj.add(new ArrayList<>());
        for (int[] e : E) { adj.get(e[0]).add(e[1]); adj.get(e[1]).add(e[0]); }

        // --- atoms in rings: shortest path between the endpoints of each closure edge ---
        boolean[] inRing = new boolean[n];
        List<Integer> ringSizes = new ArrayList<>();
        int nAromRing = 0;
        for (int[] c : CLOSE) {
            List<Integer> path = shortestPath(adj, c[0], c[1], c);
            if (path == null) continue;
            ringSizes.add(path.size());
            boolean allArom = true;
            for (int v : path) { inRing[v] = true; if (!V.get(v).isAromatic()) allArom = false; }
            if (allArom) nAromRing++;
        }

        // --- atom counts (exclude MAPStar: it is the backbone carbon, not part of the modification) ---
        int nC=0,nO=0,nN=0,nS=0,nP=0,nHal=0,nArom=0,nStar=0,nAtom=0;
        for (MAPAtomAbstract a : V) {
            if (a instanceof MAPStar) { nStar++; continue; }
            nAtom++;
            String s = a.getSymbol();
            if (a.isAromatic()) nArom++;
            switch (s) {
                case "C": nC++; break;  case "O": nO++; break;  case "N": nN++; break;
                case "S": nS++; break;  case "P": nP++; break;
                case "F": case "Cl": case "Br": case "I": nHal++; break;
            }
        }

        // --- carbonyls: C with a double bond to O ---
        int nCarbonyl = 0, nCarboxylLike = 0;
        for (int i = 0; i < n; i++) {
            if (!"C".equals(V.get(i).getSymbol())) continue;
            int dblO = 0, sglO = 0;
            for (int[] e : E) {
                int o = (e[0]==i) ? e[1] : (e[1]==i) ? e[0] : -1;
                if (o < 0) continue;
                if (!"O".equals(V.get(o).getSymbol())) continue;
                if (e[2] == 2) dblO++; else sglO++;
            }
            if (dblO >= 1) nCarbonyl++;
            if (dblO >= 1 && sglO >= 1) nCarboxylLike++;   // ester / carboxyl
        }

        // --- longest acyclic carbon chain (longest path in the non-ring carbon forest) ---
        int maxChain = longestCarbonChain(adj, inRing);

        // --- attachment atom: symbol of the star's neighbour ---
        String attach = "?";
        for (int i = 0; i < n; i++) {
            if (!(V.get(i) instanceof MAPStar)) continue;
            for (int o : adj.get(i)) { attach = V.get(o).getSymbol(); break; }
            break;
        }

        StringBuilder rs = new StringBuilder();
        Collections.sort(ringSizes);
        for (int i = 0; i < ringSizes.size(); i++) { if (i>0) rs.append(","); rs.append(ringSizes.get(i)); }

        return String.join("\t", map, rest, "OK",
            ""+nAtom, ""+nC, ""+nO, ""+nN, ""+nS, ""+nP, ""+nHal,
            ""+nArom, ""+CLOSE.size(), ""+nAromRing, rs.length()==0?"-":rs.toString(),
            ""+nCarbonyl, ""+nCarboxylLike, ""+maxChain, attach, ""+nStar);
    }

    /** register vertices and edges, resolving MAPAtomCyclic into ring-closure edges */
    static void collect(MAPGraph g) {
        for (MAPAtomAbstract a : g.getAtoms())
            if (!(a instanceof MAPAtomCyclic) && !IDX.containsKey(a)) { IDX.put(a, V.size()); V.add(a); }
        for (MAPGraph c : g.getChildGraphs()) collect(c);
        for (MAPAtomAbstract a : g.getAtoms()) {
            if (a instanceof MAPAtomCyclic) continue;
            Integer u = IDX.get(a);
            if (u == null) continue;
            for (MAPConnection cn : a.getChildConnections()) {
                MAPAtomAbstract t = cn.getAtom();
                int order = cn.getBondType() == null ? 1 : Math.max(cn.getBondType().getNumber(), 1);
                if (t instanceof MAPAtomCyclic) {
                    Integer v = IDX.get(((MAPAtomCyclic) t).getCyclicAtom());
                    if (v != null) { E.add(new int[]{u, v, order}); CLOSE.add(new int[]{u, v, order}); }
                } else {
                    Integer v = IDX.get(t);
                    if (v != null) E.add(new int[]{u, v, order});
                }
            }
        }
        for (MAPGraph c : g.getChildGraphs()) collectEdges(c);
    }
    static void collectEdges(MAPGraph g) { /* edges already handled in collect() recursion */ }

    /** shortest path from s to t without using the given closure edge */
    static List<Integer> shortestPath(List<List<Integer>> adj, int s, int t, int[] skip) {
        int n = adj.size();
        int[] prev = new int[n]; Arrays.fill(prev, -2);
        ArrayDeque<Integer> q = new ArrayDeque<>();
        q.add(s); prev[s] = -1;
        boolean skipped = false;
        while (!q.isEmpty()) {
            int u = q.poll();
            for (int v : adj.get(u)) {
                if (!skipped && ((u==skip[0]&&v==skip[1])||(u==skip[1]&&v==skip[0]))) { skipped = true; continue; }
                if (prev[v] != -2) continue;
                prev[v] = u; q.add(v);
            }
        }
        if (prev[t] == -2) return null;
        List<Integer> path = new ArrayList<>();
        for (int c = t; c != -1; c = prev[c]) path.add(c);
        return path;
    }

    /** longest path in the subgraph of acyclic carbons (a forest -> diameter per component) */
    static int longestCarbonChain(List<List<Integer>> adj, boolean[] inRing) {
        int n = adj.size();
        boolean[] ok = new boolean[n];
        for (int i = 0; i < n; i++)
            ok[i] = "C".equals(V.get(i).getSymbol()) && !inRing[i] && !(V.get(i) instanceof MAPStar);
        boolean[] seen = new boolean[n];
        int best = 0;
        for (int i = 0; i < n; i++) {
            if (!ok[i] || seen[i]) continue;
            int[] f1 = bfsFar(adj, ok, seen, i, true);
            int[] f2 = bfsFar(adj, ok, new boolean[n], f1[0], false);
            best = Math.max(best, f2[1] + 1);
        }
        return best;
    }
    static int[] bfsFar(List<List<Integer>> adj, boolean[] ok, boolean[] mark, int s, boolean record) {
        int n = adj.size();
        int[] d = new int[n]; Arrays.fill(d, -1);
        ArrayDeque<Integer> q = new ArrayDeque<>();
        q.add(s); d[s] = 0; if (record) mark[s] = true;
        int far = s;
        while (!q.isEmpty()) {
            int u = q.poll();
            if (d[u] > d[far]) far = u;
            for (int v : adj.get(u)) {
                if (!ok[v] || d[v] != -1) continue;
                d[v] = d[u] + 1; if (record) mark[v] = true; q.add(v);
            }
        }
        return new int[]{far, d[far]};
    }
}
