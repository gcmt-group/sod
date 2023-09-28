BEGIN {
  print "Space Group number",spgr ;
}


 {
  for (i=0;i<noperators;i++) 
   {
    for(j=0;j<3;j++) {
      operator[i,j]=$(j+1)}
      getline;
    }
  }


 {
  for (i=0;i<noperators;i++) 
   {
    print i+1;
    for(j=0;j<3;j++) {

     if      (operator[i,j]=="x")      {print "1   0   0  0"}  
     else if (operator[i,j]=="y")      {print "0   1   0  0"}  
     else if (operator[i,j]=="z")      {print "0   0   1  0"}  
     else if (operator[i,j]=="-x")     {print "-1  0   0  0"}  
     else if (operator[i,j]=="-y")     {print "0  -1   0  0"}  
     else if (operator[i,j]=="-z")     {print "0   0  -1  0"}  

     else if (operator[i,j]=="x+1/2")  {print "1   0   0  0.5"}  
     else if (operator[i,j]=="y+1/2")  {print "0   1   0  0.5"}  
     else if (operator[i,j]=="z+1/2")  {print "0   0   1  0.5"}  
     else if (operator[i,j]=="-x+1/2") {print "-1  0   0  0.5"}  
     else if (operator[i,j]=="-y+1/2") {print "0  -1   0  0.5"}  
     else if (operator[i,j]=="-z+1/2") {print "0   0  -1  0.5"}  

     else if (operator[i,j]=="x+1/3")  {print "1   0   0  0.33"}  
     else if (operator[i,j]=="y+1/3")  {print "0   1   0  0.33"}  
     else if (operator[i,j]=="z+1/3")  {print "0   0   1  0.33"}  
     else if (operator[i,j]=="-x+1/3") {print "-1  0   0  0.33"}  
     else if (operator[i,j]=="-y+1/3") {print "0  -1   0  0.33"}  
     else if (operator[i,j]=="-z+1/3") {print "0   0  -1  0.33"}  

     else if (operator[i,j]=="x+2/3")  {print "1   0   0  0.66"}  
     else if (operator[i,j]=="y+2/3")  {print "0   1   0  0.66"}  
     else if (operator[i,j]=="z+2/3")  {print "0   0   1  0.66"}  
     else if (operator[i,j]=="-x+2/3") {print "-1  0   0  0.66"}  
     else if (operator[i,j]=="-y+2/3") {print "0  -1   0  0.66"}  
     else if (operator[i,j]=="-z+2/3") {print "0   0  -1  0.66"}  

     else if (operator[i,j]=="x+1/4")  {print "1   0   0  0.25"}  
     else if (operator[i,j]=="y+1/4")  {print "0   1   0  0.25"}  
     else if (operator[i,j]=="z+1/4")  {print "0   0   1  0.25"}  
     else if (operator[i,j]=="-x+1/4") {print "-1  0   0  0.25"}  
     else if (operator[i,j]=="-y+1/4") {print "0  -1   0  0.25"}  
     else if (operator[i,j]=="-z+1/4") {print "0   0  -1  0.25"}  

     else if (operator[i,j]=="x+3/4")  {print "1   0   0  0.75"}  
     else if (operator[i,j]=="y+3/4")  {print "0   1   0  0.75"}  
     else if (operator[i,j]=="z+3/4")  {print "0   0   1  0.75"}  
     else if (operator[i,j]=="-x+3/4") {print "-1  0   0  0.75"}  
     else if (operator[i,j]=="-y+3/4") {print "0  -1   0  0.75"}  
     else if (operator[i,j]=="-z+3/4") {print "0   0  -1  0.75"}  

     else  print "ERROR in Group",spgr,"Operator",i+1,"Axis",j+1,operator[i,j]}
    }
  }


