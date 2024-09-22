#  awk  program to convert a string of eight decimal digits
#  to a simple quantum register expression.
#  FIT, Monash University, Sept 2024
{
  split($0, digits, "")
  n = length($0);
  if (n != 8)
  {
    printf("ID should have eight digits.\n");
  }
  else
  {
    for (i=3; i >= 1; i--)
    {
      r = digits[i] % 4;
      if  (r == 0)
      {
        factor[i] = "H";
      }
      else  if  (r == 1)
      {
        factor[i] = "X";
      }
      else  if  (r == 2)
      {
        factor[i] = "Y";
      }
      else  #  r == 3
      {
        factor[i] = "Z";
      }
    }
    s = int(digits[4] / 2) % 4;
    if  (s == 0)
    {
      aFactor = "H";
    }
    else  if  (s == 1)
    {
      aFactor = "X";
    }
    else  if  (s == 2)
    {
      aFactor = "Y";
    }
    else  #  s == 3
    {
      aFactor = "Z";
    }
    if  (digits[4] % 2 == 1)
    {
      factor[4] = "CNOT";
      factor[5] = aFactor; 
    }
    else
    {
      factor[4] = aFactor;
      factor[5] = "CNOT"; 
    }
    for (i=7; i >= 5; i--)
    {
      r = digits[i] % 4;
      if  (r == 0)
      {
        factor[i+1] = "H";
      }
      else  if  (r == 1)
      {
        factor[i+1] = "X";
      }
      else  if  (r == 2)
      {
        factor[i+1] = "Y";
      }
      else  #  r == 3
      {
        factor[i+1] = "Z";
      }
    }
    factor[9] = "k" (int(digits[8] / 4) % 2);
    factor[10] = "k" (int(digits[8] / 2) % 2);
    factor[11] = "k" (digits[8] % 2);
    printf(" (%s (x) %s (x) %s) * (%s (x) %s) * (%s (x) %s (x) %s) * (%s (x) %s (x) %s)\n", factor[1], factor[2], factor[3], factor[4], factor[5], factor[6], factor[7], factor[8], factor[9], factor[10], factor[11]);
  }
}
