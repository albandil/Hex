/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_DB_INTERFACES
#define HEX_DB_INTERFACES

extern "C"
{

/**
 * @brief Initialize the environment.
 * 
 * The function will open existing database with path specified in "dbname".
 * If the database does not exist, a new empty database will be created.
 * Also, all registered variables will be given chance to initialize by
 * the function
 * @code
 *    Variable::initialize (sqlitepp::session & db)
 * @endcode
 * This can be used to register special functions needed later in the SQL
 * statements. For example the variable IntegralCrossSection registeres
 * two functions: the simple square root and a more complicated numerical
 * integrator that uses ClenshawCurtis quadrature.
 * 
 * The function "initialize" also disables journalling to speed up insertion
 * of data.
 */
void initialize (const char* dbname);

/**
 * @brief Create new database.
 * 
 * This will create the necessary table structure in a new empty database.
 * It is necessary to call @ref initialize first, so that the databse file
 * is created and opened. The tables are created by the call to
 * @code
 *     Variable::SQL_CreateTable();
 * @endcode
 * for every registered variable. The return value is a vector of SQL statement
 * which are executed.
 */
void create_new_database ();

/**
 * @brief Import SQL batch file.
 * 
 * The effect of this function should be equivalent to direct use of sqlite3
 * program:
 * @code
 *     sqlite3 hex.db < batchfile.sql
 * @endcode
 * 
 * The function uses one line at a time and splits the lines into statements
 * on semicolons. For this reason it does not allow multi-line statements.
 * The statements produced by computational modules have always one statement
 * per line.
 */
void import (const char* sqlname);

/**
 * @brief Update integral and complete cross section.
 * 
 * This function is used after import of T-matrices. All variables will be
 * given chance to update themselves using the function
 * @code
 *     Variable::SQL_Update();
 * @endcode
 * which returns SQL statement that should do the update. If a variable
 * uses other variables besides T-matrices for its update, it must be inserted
 * to VariableList AFTER those variables, because the order of update is
 * set by the order in the list.
 */
void update ();

/**
 * @brief Optimize the SQLite database file.
 * 
 * Reduces the occupied space using the VACUUM statement.
 */
void optimize ();

/**
 * @brief Dump contents of the T-matrix table.
 * 
 * The output can be used to construct an equivalent table.
 * The corresponding code is
 * @code
 * sqlite> .mode insert
 * sqlite> select * from tmat
 * @endcode
 */
void dump (const char* dumpname);

/**
 * @brief Scattering anplitude.
 * 
 * @param ni Initial atomic principal quantum number.
 * @param li Initial atomic orbital quantum number.
 * @param mi Initial atomic magnetic quantum number.
 * @param nf Final atomic principal quantum number.
 * @param lf Final atomic orbital quantum number.
 * @param mf Final atomic magnetic quantum number.
 * @param S Total spin (0 = singlet, 1 = triplet).
 * @param E Impact energy in Rydbergs.
 * @param N Sample count.
 * @param angles Real array of length N containing scattering angles.
 * @param result Complex array of length N (or real array of length 2N) to contain the amplitudes.
 */
void scattering_amplitude
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S, double E,
    int N,
    double const * const angles,
    double  * const result
);

/**
 * @brief Differential cross section.
 * 
 * @param ni Initial atomic principal quantum number.
 * @param li Initial atomic orbital quantum number.
 * @param mi Initial atomic magnetic quantum number.
 * @param nf Final atomic principal quantum number.
 * @param lf Final atomic orbital quantum number.
 * @param mf Final atomic magnetic quantum number.
 * @param S Total spin (0 = singlet, 1 = triplet).
 * @param E Impact energy in Rydbergs.
 * @param N Sample count.
 * @param angles Real array of length N containing scattering angles.
 * @param result Real array of length N to contain the cross sections.
 */
void differential_cross_section
(
    int ni, int li, int mi,
    int nf, int lf, int mf,
    int S,
    double E,
    int N,
    double const * const angles,
    double       * const dcs
);

};

#endif
