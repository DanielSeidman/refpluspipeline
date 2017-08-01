/**
* Author: Jan Schroeder
* Modified by Marek Cmero
**/

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;


public class ReferenceIndex {

	Hashtable<String, Integer> str_to_int;
	ArrayList<String> int_to_str;
	ArrayList<Integer> chr_lengths;

	ReferenceIndex(String index_filename) throws IOException{
		str_to_int = new Hashtable<String, Integer>();
		int_to_str = new ArrayList<String>();
		chr_lengths = new ArrayList<Integer>();

		BufferedReader in = new BufferedReader(new FileReader(index_filename));
		String line;
		int i = 0;
		while ( (line = in.readLine()) != null) {
			StringTokenizer t = new StringTokenizer(line);
			String name = t.nextToken();
			int len = Integer.parseInt(t.nextToken());
			str_to_int.put(name, i);
			int_to_str.add(name);
			chr_lengths.add(len);

			i++;
		}
	}

	public int getLength(String name){
		return getLength(str_to_int.get(name));
	}
	public int getLength(int index){
		return chr_lengths.get(index);
	}
	public int getIndex(String name){
		return str_to_int.get(name);
	}
	public String getChromosomeName(int index){
		return int_to_str.get(index);
	}
	public int size(){
		return int_to_str.size();
	}
}
