# A simple way to patch qrefine into PHENIX
git clone https://github.com/qrefine/qr-core.git
mv qr-core core
phenix.python -m pip install -r requirements.txt
phenix.python -m pip install -r core/requirements.txt
cd core/plugin
git clone https://github.com/qrefine/qr-plugin-yoink.git
git clone https://github.com/qrefine/qr-plugin-ase.git
cd ..
cd tests
git clone https://github.com/qrefine/qr-tests-p1.git
git clone https://github.com/qrefine/qr-tests-cluster.git
libtbx.configure qrefine
