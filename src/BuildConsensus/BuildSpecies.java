package BuildConsensus;

import jaligner.Sequence;
import jaligner.matrix.MatrixLoaderException;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;


public class BuildSpecies {

	private static Sequence seqLoaded;
	private static List<SAMRecord> records;
	private static BuildClasses holder;
	private static HashMap<Integer,BuildClasses.InDelInfo> deletions;
	private static HashMap<Integer,BuildClasses.InDelInfo> insertions;
	private static List<BuildClasses.cutter> toCut;
	private static int start;
	private static int end;
	
	public static String getSpecies (Sequence seqLoaded, List<SAMRecord> records, int start, int end) throws IOException, InterruptedException, MatrixLoaderException, URISyntaxException {
		
		BuildSpecies.seqLoaded = seqLoaded;
		BuildSpecies.records = records;
		BuildSpecies.start = start;
		BuildSpecies.end = end;
		holder = new BuildClasses();
		deletions = new HashMap<Integer,BuildClasses.InDelInfo>();
		insertions = new HashMap<Integer,BuildClasses.InDelInfo>();	
		toCut  = new ArrayList<BuildClasses.cutter>();

		String consensus = startSpecies();
		
		return consensus;
		
	}
	
	private static String startSpecies() throws IOException, InterruptedException, MatrixLoaderException, URISyntaxException {
		
		BuildQuals build = new BuildQuals();
		HashMap<Integer,String> mutations = new LinkedHashMap<Integer,String>();

		//Populate the map, so that I can add stuff to it
		for(int parser = start; parser <= end; parser++) {
			mutations.put(parser, "");
		}
		
		//Now iterate over the readblocks for each read!
		for (SAMRecord thisOne : records) {
			
			String readSequence = thisOne.getReadString();
			char[] seqSubArray = readSequence.toCharArray();
			char[] seqQualArray = thisOne.getBaseQualityString().toCharArray();
			
			List<BuildClasses.AlignmentBlob> blobs = getReadBlocks(thisOne.getCigar(), thisOne.getAlignmentStart(), thisOne.getReadName());			
			
			for (BuildClasses.AlignmentBlob blo : blobs) {
				
				for (int refParser = blo.getReferenceStart(), seqParser = blo.getReadStart() - 1; refParser <= ((blo.getReferenceStart() + blo.getLength()) - 1); refParser++, seqParser++) {
					if (refParser < start || refParser > end) {
						continue;
					}
					if (build.qualityString(String.valueOf(seqQualArray[seqParser])) > 10.0) {
						String newString = mutations.get(refParser) + seqSubArray[seqParser];
						newString = newString.toUpperCase();
						mutations.put(refParser, newString);
					}

				}
				
			}
		
		}
		
		// Generate a pileup of all the read information at this site, and decide consensus (all the while caring about ins/del status):
		String consensus = getPileup(mutations);
		
		// Check for deletions in the consensus:
		checkDel(consensus);
			
		// Check for insertions in the consensus:
		checkIns(consensus);
		
		// Now apply InDels to the consensus:
		// Have to apply the change to the consensus, and then subtract or add to each element in Array
		
		String meiCut = seqLoaded.getSequence();
		for (int x = 0; x < toCut.size(); x++) {
			BuildClasses.cutter current = toCut.get(x);
			if (current.getType().equals("insertion")) {
				String before = consensus.substring(0, current.getLocation());
				String after = consensus.substring(current.getLocation());
				int length = current.getSequence().length();
				consensus = before + current.getSequence() + after;
				before = meiCut.substring(0, current.getLocation());
				after = meiCut.substring(current.getLocation());
				length = current.getSequence().length();
				meiCut = before + current.getSequence() + after;
				// Now go through anything after and add length to every position
				for (int y = x+1; y < toCut.size(); y++) {
					BuildClasses.cutter fixed = toCut.get(y);
					fixed.setLocation(fixed.getLocation() + length);
				}
			} else {
				String before = consensus.substring(0, current.getLocation());
				String after = consensus.substring(current.getLocation() + current.getSequence().length(), meiCut.length());
				int length = current.getSequence().length();
				consensus = before + after;
				before = meiCut.substring(0, current.getLocation());
				after = meiCut.substring(current.getLocation() + current.getSequence().length(), meiCut.length());
				length = current.getSequence().length();
				meiCut = before + after;
				// Now go through anything after and subtract length from every position
				for (int y = x+1; y < toCut.size(); y++) {
					BuildClasses.cutter fixed = toCut.get(y);
					fixed.setLocation(fixed.getLocation() - length);
				}
			}
		}
		
		return consensus;
		
	}

	private static List<BuildClasses.AlignmentBlob> getReadBlocks(Cigar cigar, int readStart, String readName) {
		
		List<BuildClasses.AlignmentBlob> blobs = new ArrayList<BuildClasses.AlignmentBlob>();
				
		int currentRefPos = readStart;
		int currentReadPos = 1;
		
		for (CigarElement ele : cigar.getCigarElements()) {
			
			if (ele.getOperator().toString().equals("S")) {
				currentReadPos += ele.getLength();
				continue;
			} else if (ele.getOperator().toString().equals("M")) {
				BuildClasses.AlignmentBlob blob = holder.new AlignmentBlob();
				blob.setLength(ele.getLength());
				blob.setReadStart(currentReadPos);
				blob.setReferenceStart(currentRefPos);
				blobs.add(blob);
				currentReadPos+= ele.getLength();
				currentRefPos+= ele.getLength();
			} else if (ele.getOperator().toString().equals("D")) {
				if (deletions.containsKey(currentRefPos)) {
					BuildClasses.InDelInfo delInfo = deletions.get(currentRefPos);
					delInfo.setTotal(delInfo.getTotal() + 1);
					deletions.put(currentRefPos, delInfo);
				} else {
					BuildClasses.InDelInfo delInfo = holder.new InDelInfo();
					delInfo.setTotal(1);
					delInfo.setLength(ele.getLength());
					deletions.put(currentRefPos, delInfo);
				}
				currentRefPos+= ele.getLength();
				continue;
			} else if (ele.getOperator().toString().equals("I")) {
				if (insertions.containsKey(currentRefPos)) {
					BuildClasses.InDelInfo insInfo = insertions.get(currentRefPos);
					insInfo.setTotal(insInfo.getTotal() + 1);
					BuildClasses.InDelInfo.Info info = insInfo.new Info();
					info.setInsLen(ele.getLength());
					info.setPosition(currentReadPos);
					insInfo.addToInsPos(readName, info);
					insertions.put(currentRefPos, insInfo);
				} else {
					BuildClasses.InDelInfo insInfo = holder.new InDelInfo();
					insInfo.setTotal(1);
					BuildClasses.InDelInfo.Info info = insInfo.new Info();
					info.setInsLen(ele.getLength());
					info.setPosition(currentReadPos);
					insInfo.addToInsPos(readName, info);
					insertions.put(currentRefPos, insInfo);
				}
				currentReadPos+=ele.getLength();
				continue;
			} else {
				continue;
			}
		}
		
		return blobs;
		
	}
	private static String getPileup (HashMap<Integer,String> mutations) {
		
		String consensus = "";
		for (int x = start; x <= end; x++) {
			
			if (mutations.get(x).length() > 0) {
								
				double aTotal = 0;
				double tTotal = 0;
				double cTotal = 0;
				double gTotal = 0;
				char bases[] =  mutations.get(x).toCharArray();
				double numBases = bases.length;
				
				for (char base : bases) {
						
					if (base == 'A') {
						aTotal++;
					} else if (base == 'T') {
						tTotal++;
					} else if (base == 'G') {
						gTotal++;
					} else {
						cTotal++;
					}
						
				}
					
				aTotal /= numBases;
				tTotal /= numBases;
				cTotal /= numBases;
				gTotal /= numBases;
					
				if (aTotal >= tTotal && aTotal >= cTotal && aTotal >= gTotal) {
					consensus += "A";
				} else if (tTotal >= aTotal && tTotal >= cTotal && tTotal >= gTotal) {
					consensus += "T";
				} else if (cTotal >= aTotal && cTotal >= tTotal && cTotal >= gTotal) {
					consensus += "C";
				} else {
					consensus += "G";
				}
				
			} else {
				consensus += "N";	
			}
			
		}
		
		return consensus;
		
	}
	private static void checkDel(String consensus) {
		
		for (Map.Entry<Integer, BuildClasses.InDelInfo> entry : deletions.entrySet()) {
			
			double totalHits = entry.getValue().getTotal();
			int delLoc = entry.getKey();
			double totalReads = 0;
			
			for (SAMRecord thisOne : records) {
				
				if (thisOne.getAlignmentStart() < delLoc && thisOne.getAlignmentEnd() > delLoc) {
					
					totalReads++;
					
				}
				
			}
			
			delLoc = delLoc - start;
			
			double calc = totalHits / totalReads;
			
			if ((calc > 0.50) && totalReads > 2) {
				
				String del = consensus.substring(delLoc, delLoc + entry.getValue().getLength());
				BuildClasses.cutter cut = holder.new cutter();
				cut.setLocation(delLoc);
				cut.setSequence(del);
				cut.setType("deletion");
				// add items in order to the array
				if (toCut.isEmpty()) {
					toCut.add(cut);	
				} else {
					int size = toCut.size();
					if (size == 1) {
						if (toCut.get(0).getLocation() > delLoc) {
							toCut.add(0, cut);
						} else {
							toCut.add(cut);
						}
					} else {
						// Walk through each array element to see where the appropriate place is
						int cutSize = toCut.size();
						for (int x = 0; x < cutSize; x++) {
							if (toCut.get(x).getLocation() > delLoc) {
								toCut.add(x, cut);
								break;
							}
							if (x == toCut.size() - 1) {
								toCut.add(cut);
							}
						}
					}
				}	
			}
		}
	}
	private static void checkIns(String consensus) {
		
		for (Map.Entry<Integer, BuildClasses.InDelInfo> entry : insertions.entrySet()) {
			
			double totalHits = entry.getValue().getTotal();
			int insLoc = entry.getKey();
			double totalReads = 0;
			Map<String, Integer> insConsensus = new HashMap<String, Integer>();
			
			for (SAMRecord thisOne : records) {
				
				if (thisOne.getAlignmentStart() < insLoc && thisOne.getAlignmentEnd() > insLoc) {
					
					totalReads++;
					if (entry.getValue().getInsPos().containsKey(thisOne.getReadName())) {
						int insReadLoc = entry.getValue().getInsPos().get(thisOne.getReadName()).getPosition();
						int insLength = entry.getValue().getInsPos().get(thisOne.getReadName()).getInsLen();
						String insString = thisOne.getReadString().substring(insReadLoc, insReadLoc + insLength);
						if (insConsensus.containsKey(insString)) {
							int adjust = insConsensus.get(insString);
							adjust++;
							insConsensus.put(insString, adjust);
						} else {
							insConsensus.put(insString, 1);
						}
					}
					
				}
				
			}
			
			double calc = totalHits / totalReads;
			
			insLoc = insLoc - start;
			
			if ((calc >= 0.50) && totalReads > 2) {
			
				//  Decide the proper insertion to place into the consensus:
				String highest = "";
				int currentScore = 0;
				
				for (Map.Entry<String, Integer> comp : insConsensus.entrySet()) {
					
					if (comp.getValue() > currentScore) {
						
						currentScore = comp.getValue();
						highest = comp.getKey();
						
					}
					
				}
				
				BuildClasses.cutter cut = holder.new cutter();
				cut.setLocation(insLoc);
				cut.setSequence(highest);
				cut.setType("insertion");
				if (toCut.isEmpty()) {
					toCut.add(cut);	
				} else {
					int size = toCut.size();
					if (size == 1) {
						if (toCut.get(0).getLocation() > insLoc) {
							toCut.add(0, cut);
						} else {
							toCut.add(cut);
						}
					} else {
						// Walk through each array element to see where the appropriate place is
						int cutSize = toCut.size();
						for (int x = 0; x < cutSize; x++) {
							
							if (toCut.get(x).getLocation() > insLoc) {
								toCut.add(x, cut);
								break;
							}
							if (x == (toCut.size() - 1)) {
								toCut.add(cut);
							}
						}
					}
				}	
			}
		}		
	}
}
