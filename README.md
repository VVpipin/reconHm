# reconHm
reconstruction vector potential  and  magnetic helicity maps 
The scripts in this repository design to reconstruct magnetic vector-potential components from
synoptic maps of the magnetic field vector maps, which are currently produced by SDO/HMI and SOLIS/VSM.
In result we can produce the magnetic helicity density maps and  B_r (curl B)_r as well. Both quantities
are rather important for the solar dynamo diagnostic. The algorithm stems from the basic decompositon
of the vector field on sum of curls. It is employed in dynamo models (see Krause & Raedler 1980).
The method can be applied only for the global decomposition of the vector magnetic field. The local analysis of this type
meets the problem of gauge.
There is a test field distributions obtained in the model presented recently  by Pipin & Kosovichev 
doi:10.3847/1538-4357/aae1fb. It is in test-nxpgln.py (it uses data stored from the dynamo model to npz files).
Note that magnetic field in that model is in dimensionless form.
I have used nxpgln.py (with auxilary routine initrmy.py) for reconstruction SDO/HMI with spherical harmonics ell=179
To get the smoothed helicity maps I have to smooth the original data as well (see inside  nxpgln.py, window size=[4,8]).
This is still very high resolution for the current helicity calculations.  I found that size=[10,20] may give
reasonable agreement between l-spectra of the current and magnetic helicity, see figure with title('croscorr AB vs hc')
in every case. The file pglnhmi360720.py uses auxilary routine initr.py and ell=299. Hope it is enough to resolve HMI structure.
Not sure about current helicity which may need the higher ell.

Note that to run these files you need the proper python distribution with axulary packages numpy, scipy, astropy(reading fits)
and PySHTools. All this available on repository of manjaro linux (derivative of archlinux). I don't know about others.
Also in test files there are some path settings to save figures. Be aware!
