#!/bin/bash

# installs fnm (Fast Node Manager)
curl -fsSL https://fnm.vercel.app/install | bash
source ~/.bashrc
# download and install Node.js
fnm use --install-if-missing 20
# verifies the right Node.js version is in the environment
node -v # should print `v20.13.1`
# verifies the right NPM version is in the environment
npm -v # should print `10.5.2`

npm install -g http-server