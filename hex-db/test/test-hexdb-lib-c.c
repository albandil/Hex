#include <stdio.h>
#include <stdlib.h>

#include "../src/interfaces.h"

int main (void)
{
    // setup the database (create new if non-existant)
    hex_initialize("hex-c.db");
    
    // create new content in the database
    hex_new();
    
    // create some data-file
    FILE* sql = fopen("batch-c.sql", "w");
    fprintf(sql, "BEGIN TRANSACTION;\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,0,0.65,0,1.370592e+01,1.210727e+01,0,0);\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,1,0.65,0,-1.050670e+00,2.758337e+01,0,0);\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,0,0.75,0,1.231381e+01,9.205694e+00,0,0);\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,1,0.75,0,8.642624e-01,2.568767e+01,0,0);\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,0,0.85,0,1.096124e+01,9.684800e+00,0,0);\n");
    fprintf(sql, "INSERT OR REPLACE INTO \"tmat\" VALUES (1,0,0,1,0,0,0,1,0.85,0,2.304937e+00,2.392586e+01,0,0);\n");
    fprintf(sql, "COMMIT;\n");
    fclose(sql);
    
    // import the data-file
    hex_import("batch-c.sql");
    
    // precompute cross sections
    hex_update();
    
    // compute the cross sections
    double energy[3] = { 0.65, 0.75, 0.85 };
    double ccs[3];
    hex_complete_cross_section
    (
        1,0,0,      // initial state (1s)
        1,0,0,      // final state (1s)
        3,          // number of energies
        energy,     // energy list
        ccs,        // cross sections
        NULL        // auxiliary pointer
    );
    
    // print the cross section
    printf("  energy                     cross section\n");
    printf("  %.17f        %.15f\n", energy[0], ccs[0]);
    printf("  %.17f        %.15f\n", energy[1], ccs[1]);
    printf("  %.17f        %.15f\n", energy[2], ccs[2]);
    
    // done
    return EXIT_SUCCESS;
}
