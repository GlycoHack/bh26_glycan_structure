import java.io.*;
import java.util.*;
import org.glycoinfo.WURCSFramework.wurcs.graph.CarbonDescriptor;

/** Structural features of a WURCS SkeletonCode, using the framework's own
 *  CarbonDescriptor table. One character = one backbone carbon; the first and
 *  last characters are terminal carbons, the rest non-terminal, and several
 *  characters mean different things in the two positions.
 *  Input: one SkeletonCode per line. Output: TSV. */
public class SkeletonFeatures {

    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));
        out.println(String.join("\t", "skeleton","status","nC","first","last",
            "nAnomer","nStereoDef","nStereoUnk","nDeoxy","nUndef","nCarbonyl","nAcid",
            "nDouble","nTriple","nSP3","nSP2","nSP","nSPX","unknownChars"));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) continue;
            String sk = line.split("\t")[0];
            out.println(features(sk));
        }
        out.flush();
    }

    static String features(String sk) {
        int n = sk.length();
        int nAnomer=0, nStereoDef=0, nStereoUnk=0, nDeoxy=0, nUndef=0;
        int nCarbonyl=0, nAcid=0, nDouble=0, nTriple=0;
        int nSP3=0, nSP2=0, nSP=0, nSPX=0;
        StringBuilder bad = new StringBuilder();
        String first = "-", last = "-";

        for (int i = 0; i < n; i++) {
            char c = sk.charAt(i);
            boolean terminal = (i == 0 || i == n - 1);
            CarbonDescriptor cd = CarbonDescriptor.forCharacter(c, terminal);
            if (cd == null) { bad.append(c); continue; }
            if (i == 0)     first = cd.name();
            if (i == n - 1) last  = cd.name();

            String orb = cd.getHybridOrbital();
            if      ("sp3".equals(orb)) nSP3++;
            else if ("sp2".equals(orb)) nSP2++;
            else if ("sp".equals(orb))  nSP++;
            else                        nSPX++;

            int b1 = cd.getBondTypeCarbon1(), b2 = cd.getBondTypeCarbon2();
            if (b1 == 3 || b2 == 3) nTriple++;
            else if (b1 == 2 || b2 == 2) nDouble++;

            String nm = cd.name();
            if (nm.contains("ANOMER"))            nAnomer++;
            if (nm.contains("STEREO") || nm.contains("CHIRAL")) {
                if (nm.endsWith("_X_L") || nm.endsWith("_X_U")) nStereoUnk++; else nStereoDef++;
            }
            if (nm.contains("METHYNE") || nm.contains("METHYL")) nDeoxy++;
            if (nm.contains("UNDEF") || nm.contains("UNKNOWN"))  nUndef++;
            if (nm.contains("ALDEHYDE") || nm.contains("KETONE"))nCarbonyl++;
            if (nm.contains("ACID"))                             nAcid++;
        }
        return String.join("\t", sk, bad.length()==0 ? "OK" : "UNKNOWN_CHAR",
            ""+n, first, last, ""+nAnomer, ""+nStereoDef, ""+nStereoUnk, ""+nDeoxy,
            ""+nUndef, ""+nCarbonyl, ""+nAcid, ""+nDouble, ""+nTriple,
            ""+nSP3, ""+nSP2, ""+nSP, ""+nSPX, bad.length()==0 ? "-" : bad.toString());
    }
}
