# A simple way to patch qrefine into PHENIX, can change to python later.
git clone https://github.com/qrefine/qr-core.git
mv qr-core core
libtbx.python -m pip install -r requirements.txt
libtbx.python -m pip install -r core/requirements.txt
cd core/plugin
git clone https://github.com/qrefine/qr-plugin-yoink.git
git clone https://github.com/qrefine/qr-plugin-ase.git
cd ..
#cd tests
#git clone https://github.com/qrefine/qr-tests-p1.git  -> moved to regression
#git clone https://github.com/qrefine/qr-tests-cluster.git -> moved to regression
libtbx.configure qrefine
