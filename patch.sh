# A simple way to patch qrefine into PHENIX
git clone https://github.com/qrefine/qr-core.git
mv qr-core core
phenix.python -m pip install -r requirements.txt
phenix.python -m pip install -r core/requirements.txt
libtbx.configure qrefine
