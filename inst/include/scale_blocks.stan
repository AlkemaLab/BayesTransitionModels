vector scale_blocks(vector raw, vector sigma, array[] int start, array[] int end) {
  int s = rows(sigma);
  vector[rows(raw)] result;

  for(i in 1:s) {
    result[start[i]:end[i]] = raw[start[i]:end[i]] * sigma[i];
  }
  return result;
}
