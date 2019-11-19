# correct Illumina Infinium probe mappings and truncate regions to 1 bp (BED 9 format)
while (<>) {
    @bed = split;
    if ($bed[5] =~ /\+/) {
        $bed[1] -= 1;
    }
    $bed[2] = $bed[1] + 1;
    $bed[6] = $bed[1];
    $bed[7] = $bed[2];
    print join("\t", @bed) . "\n";
}
