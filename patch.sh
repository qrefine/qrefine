# A simple way to patch qrefine into PHENIX, can change to python later.
git clone https://github.com/qrefine/qr-core.git
mv qr-core core
libtbx.python -m pip install -r requirements.txt
libtbx.python -m pip install -r core/requirements.txt
cd core/plugin
git clone https://github.com/qrefine/qr-plugin-yoink.git
mv qr-plugin-yoink yoink
git clone https://github.com/qrefine/qr-plugin-ase.git
mv qr-plugin-ase ase
cd ../../
libtbx.configure qrefine
