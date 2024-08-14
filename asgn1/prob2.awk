BEGIN {
    first = 1;
}
# The next line of code gets the number of vertices from the first line
# and stores it in the variable n
NR == 1 { 
    n = $1; 
}

# Generate clauses ensuring that each edge must be included in the vertex cover
NR > 1 {
    if (!first) {
        printf " & ";
    }
    printf "(v%s|v%s)", $1, $3;

    # Mark this as not the first clause anymore
    first = 0;

    # Store the edge in an associative array for checking later
    a[$1, $3] = 1;
    a[$3, $1] = 1;
}

END {
    # Looping through the triples of vertices in the array we created
    for (i = 1; i <= n-2; i++) {
        for (j = i+1; j <= n-1; j++) {
            if (a[i, j] == 1) { # edge found, now check for triangles
                for (k = j+1; k <= n; k++) {
                    if (a[i, k] == 1 && a[j, k] == 1) {
                        # Triangle found! 
                        printf " & (~v%s|~v%s|~v%s)", i, j, k;
                    }
                }
            }
        }
    }
}
