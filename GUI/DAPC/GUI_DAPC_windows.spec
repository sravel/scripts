
# -*- mode: python -*-

block_cipher = None

added_files = [ ('E:\\Bayer\\ScriptsSEB\\scripts\\GUI\\DAPC\\includes\\*', 'includes'), ]

a = Analysis(['GUI_DAPC.py'],
             pathex=['E:\\Bayer\\ScriptsSEB\\scripts\\GUI\\DAPC\\'],
             binaries=None,
             datas=added_files,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='GUI_DAPC_windows',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='E:\Bayer\ScriptsSEB\scripts\GUI\DAPC\includes\icon.ico', resources=['E:\\Bayer\\ScriptsSEB\\scripts\\GUI\\DAPC\\includes\\dapc.ui'])


