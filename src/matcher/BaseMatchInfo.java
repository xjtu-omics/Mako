/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matcher;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author jiadonglin
 */
public class BaseMatchInfo {
        String  strIdentity;
        String curMatchStr;
//        List<Integer> matchedStartAt = new ArrayList<>();
        int matchedStartAt = -1;
        int baseIdx;
        char base;
        boolean isForward;
        boolean isBoth;
        List<StringIdentity> baseStringIdentitys = new ArrayList<>();
        
        public BaseMatchInfo(String base){
            this.curMatchStr = base;
        }
        
        public BaseMatchInfo(String matchStr, int matchedAt){
            this.curMatchStr = matchStr;
            this.matchedStartAt = matchedAt;
        }
//        public BaseMatchInfo(String strID, String matchStr, List<Integer> matchedAt){
//            this.curMatchStr = matchStr;
//            this.matchedStartAt.addAll(matchedAt);
//            this.strIdentity = strID;
//        }
        public BaseMatchInfo(char base, boolean isForward, boolean isBoth){
            this.base = base;
            this.isForward = isForward;
            this.isBoth = isBoth;
        }
        
        public void setMatchedAt(int idx){
            matchedStartAt = idx;
        }
        public void setBase(char base){
            this.base = base;
        }

        public void setIsBoth(boolean isBoth){
            this.isBoth = isBoth;
        }
        public void setIsForward(boolean isForward){
            this.isForward = isForward;
        }
        public void setBaseStringIDs(List<StringIdentity> identitys){
            baseStringIdentitys = identitys;
        }
        public BaseMatchInfo append(BaseMatchInfo info){
            
            matchedStartAt = info.matchedStartAt;
            String newStr = curMatchStr + info.curMatchStr;
            return new BaseMatchInfo(newStr, matchedStartAt);
        }
        public int curMatchedStrLen(){
            return curMatchStr.length();
        }
}
