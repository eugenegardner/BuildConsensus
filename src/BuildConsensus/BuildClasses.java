package BuildConsensus;

import java.util.HashMap;

public class BuildClasses {

	public class InDelInfo {
		
		private double total;
		private int length;
		private HashMap<String, BuildClasses.InDelInfo.Info> insPos;
		
		public InDelInfo() {
			
			insPos = new HashMap<String, BuildClasses.InDelInfo.Info>();
			
		}
		
		public class Info {
			
			private int InsLen;
			private int Position;
			
			public Info () {
				super();				
			}

			public int getInsLen() {
				return InsLen;
			}

			public void setInsLen(int insLen) {
				InsLen = insLen;
			}

			public int getPosition() {
				return Position;
			}

			public void setPosition(int position) {
				Position = position;
			}
			
		}
		
		public HashMap<String, BuildClasses.InDelInfo.Info> getInsPos() {
			return insPos;
		}
		public void addToInsPos(String readName, BuildClasses.InDelInfo.Info info) {
			this.insPos.put(readName, info);
		}
		public double getTotal() {
			return total;
		}
		public void setTotal(double total) {
			this.total = total;
		}
		public int getLength() {
			return length;
		}
		public void setLength(int length) {
			this.length = length;
		}
		
	}
	public class cutter {
	
		private String type;
		private String sequence;
		private int location;
		
		public String getType() {
			return type;
		}
		public void setType(String type) {
			this.type = type;
		}
		public String getSequence() {
			return sequence;
		}
		public void setSequence(String sequence) {
			this.sequence = sequence;
		}
		public int getLocation() {
			return location;
		}
		public void setLocation(int location) {
			this.location = location;
		}
		
	}
	public class AlignmentBlob {
		
		private int Length;
		private int ReadStart;
		private int ReferenceStart;
		
		public int getLength() {
			return Length;
		}
		public void setLength(int length) {
			Length = length;
		}
		public int getReadStart() {
			return ReadStart;
		}
		public void setReadStart(int readStart) {
			ReadStart = readStart;
		}
		public int getReferenceStart() {
			return ReferenceStart;
		}
		public void setReferenceStart(int referenceStart) {
			ReferenceStart = referenceStart;
		}
		
		
	}

}
