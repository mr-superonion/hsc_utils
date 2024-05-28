# Setups for camb.hmcode

Note, the default setup of pyhmcode is different from the setup in camb. We
need to use [this
setup](https://github.com/mr-superonion/utils_cosmosis/blob/b50e756c22f863f77316b844c116407fc3145cfe/pyhalofit/halofit_interface.py#L15)
to be consistent (within 1%) with the hmcode implemented in camb. Otherwise the relative bias can reach to 4% at k>30 h/Mpc for extremely small $\Omega_M\sim0.15$.
See the figures below for comparisons with camb.hmcode for the default pyhmcode setup (left) and the reconfigured pyhmcode (right):

<img src="https://user-images.githubusercontent.com/12228372/204708535-8fdb3058-b576-468a-9bd6-94a8adb58757.png" width="400">  <img src="https://user-images.githubusercontent.com/12228372/204699704-122e4503-9a21-4fb0-9d23-f6a7d2b53d82.png" width="400">
