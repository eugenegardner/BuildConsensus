package BuildConsensus;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;

public class BuildQuals {
	
	private Hashtable<Character, Integer> qualityTable;
	
	public BuildQuals() throws IOException {
		
		BufferedReader qualityStream;
		qualityTable = new Hashtable<Character, Integer>();
		
		qualityStream = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/ImportFiles/qualities.txt"), "UTF-8"));

		String myLine;
		String data[];
		
		try {
			while((myLine = qualityStream.readLine()) != null) {
				data = myLine.split("\t");
				qualityTable.put(data[0].charAt(0), Integer.parseInt(data[1]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		qualityStream.close();
		
	}
	
	public double qualityString(String errorString) {
		
		double quality = 0;
		double totalQuals = 0;
		char charArray[] = errorString.toCharArray();
		
		for (char itr : charArray) {
			totalQuals++;
			quality+=qualityTable.get(itr);
		}

		quality/=totalQuals;
		return quality;
	}
	
}
