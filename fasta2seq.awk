#################################################
# convert fasta to seq (one line)
#################################################
{
  if(substr($1,1,1)==">"&&FNR>1)
    printf("\n");
  else
    printf("%s",$1)p;
}
END{
  printf("\n");
}