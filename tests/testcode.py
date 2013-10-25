from eos import all
from tov import all
import matplotlib.pyplot as plt


TestEOS = SpecifyCore(13.7,2.7,3.0,2.0)

# what does a star look like with central energy density of 5e14 ?
TestProfile = profile(TestEOS, 5e14, 40)
energydensity = TestProfile[0]
radius = TestProfile[1]
enclosedmass = TestProfile[2] *Solarmass_km # convert units to normal ones

plt.figure()
plt.plot(radius, energydensity)
plt.xlabel('radius (km)')
plt.ylabel('energy density (g/cm^3)')
plt.title('Star with central energy 5e14 g/cm^3')
plt.show()

plt.figure()
plt.semilogy(radius, energydensity)
plt.xlabel('radius (km)')
plt.ylabel('energy density (g/cm^3)')
plt.title('Star with central energy 5e14 g/cm^3')
plt.show()


plt.plot(radius, enclosedmass)
plt.xlabel('radius (km)')
plt.ylabel('enclosed mass (Solar masses)')
plt.title('Star with central energy 5e14 g/cm^3')
plt.show()

print(' Central energy 5e14 means')
print('mass = ' + str(enclosedmass[-1]) + ' solar masses')
print('radius = ' + str(radius[-1]) + ' km')

TestMassRadius = mrarray(TestEOS, 1.5e14, 2.5e15, 50)
#confusingly, I already converted to solar masses for this function. sorry
# now the radius here is the outer radius over many star models
masses = TestMassRadius[:,2]
outerradii = TestMassRadius[:,1]
plt.plot(masses,outerradii)
plt.xlabel('Outer Radius (km)')
plt.ylabel('Total mass of star(Solar masses)')
plt.title('Properties of stars with this EOS')
plt.show()
