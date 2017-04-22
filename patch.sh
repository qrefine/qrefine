# A simple way to patch qrefine into PHENIX.
# Note: we can change this patch.sh file to a python solution later.
git clone https://github.com/qrefine/qr-core.git
mv qr-core core
libtbx.python -m pip install -r requirements.txt
libtbx.python -m pip install -r core/requirements.txt
cd core/plugin
git clone https://github.com/qrefine/qr-plugin-yoink.git
mv qr-plugin-yoink yoink
git clone https://github.com/qrefine/qr-plugin-ase.git
mv qr-plugin-ase ase
git clone https://github.com/qrefine/community
cd ../../
libtbx.configure qrefine
