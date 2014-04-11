package Recursion;

import java.util.ArrayList;
import java.util.Arrays;

public class CombinationSumII {
	
	static ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
    
    public ArrayList<ArrayList<Integer>> combinationSum2(int[] num, int target) {
        if(num.length > 0){
        	Arrays.sort(num);
            ArrayList<Integer> curRls =  new ArrayList<Integer>();
            findAllCombinations(num, target, 0, curRls);
        }
        return results;
    }
    
    public void findAllCombinations(int[] num, int remainTarget, int start, ArrayList<Integer> cur){
        if(remainTarget < 0){
            return;
        }
        else if(remainTarget == 0){
            results.add((ArrayList<Integer>)cur.clone());
            return;
        }
        else{
            for(int i = start; i < num.length; i++){
                if(i > start && num[i] == num[i - 1]){
                    continue;
                }
                cur.add(num[i]);
                findAllCombinations(num, remainTarget - num[i], i + 1, cur);
                cur.remove(cur.size() - 1);
            }
        }
    }
    
    public static void main(String args[]){
    	CombinationSumII csii = new CombinationSumII();
    	csii.combinationSum2(new int[]{1}, 1);
    	System.out.println(results);
    }
}
