//
// Created by a on 2019/10/25.
//

ofstream outdata;
outdata.precision(6);
outdata.open("test.dat", ios::out);
outdata << "Crystalline Cu atoms\n\n";
outdata << atoms_position.size() << " atoms\n";
outdata << seed_number << " atom types\n";
outdata << fixed;
outdata << "0" << " " << "500" << " xlo xhi\n";
outdata << "0" << " " << "500" << " ylo yhi\n";
outdata << "0" << " " << "500" << " zlo zhi\n";
outdata << "\n";
outdata << "Atoms\n\n";
for (int number_atom = 0; number_atom < atoms_position.size(); ++number_atom) {
outdata << number_atom + 1 << " " << (int) atoms_position[number_atom][0];
for (int number_out = 1; number_out < 4; ++number_out) {
outdata << " " << atoms_position[number_atom][number_out];
}
outdata << endl;
}
outdata.close();
cout << endl;