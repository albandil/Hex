program test
integer :: NULL
integer, parameter :: outunit=44
double precision, dimension(3) :: energy, sigma

! setup the database (create new if non-existant)
call hex_initialize('hex-f.db'//achar(0)) ! add terminating zero

! create new content in the database
call hex_new()

! create some data-file
open(unit=outunit, file='batch-f.sql', action='write', status='replace')
write(outunit,*) "BEGIN TRANSACTION;"
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,0,0.65,0,1.370592e+01,1.210727e+01,0,0);'
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,1,0.65,0,-1.050670e+00,2.758337e+01,0,0);'
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,0,0.75,0,1.231381e+01,9.205694e+00,0,0);'
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,1,0.75,0,8.642624e-01,2.568767e+01,0,0);'
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,0,0.85,0,1.096124e+01,9.684800e+00,0,0);'
write(outunit,*) 'INSERT OR REPLACE INTO "tmat" VALUES (1,0,0,1,0,0,0,1,0.85,0,2.304937e+00,2.392586e+01,0,0);'
write(outunit,*) 'COMMIT;'
close(outunit)

! import the data-file
call hex_import('batch-f.sql'//achar(0)) ! add terminating zero

! precompute cross sections
call hex_update()

! compute the cross sections
energy(1) = 0.65d+0
energy(2) = 0.75d+0
energy(3) = 0.85d+0
call hex_complete_cross_section&
(&
    1,0,0,      & ! initial state (1s)
    1,0,0,      & ! final state (1s)
    3,          & ! number of energies
    energy,     & ! energy list
    sigma,      & ! cross sections
    NULL        & ! auxiliary variable
)

! print the cross section
print*, ' energy                     cross section'
print*, energy(1), sigma(1)
print*, energy(2), sigma(2)
print*, energy(3), sigma(3)

end
