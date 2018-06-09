conv(record, l) =
{
   [n, k, i, d, aplist, charlist] = record;
   print([l, n, k, aplist, charlist]);
}

convert(records, l) = apply(r -> conv(r, l), records);
