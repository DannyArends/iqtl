import std.stdio;
import std.file;
import core.vararg;
import std.c.stdlib;

import std.c.windows.windows;
import core.sys.windows.dll;

__gshared HINSTANCE g_hInst;

import rlib;

extern (Windows){
  BOOL DllMain(HINSTANCE hInstance, ULONG ulReason, LPVOID pvReserved){
    switch (ulReason){
    case DLL_PROCESS_ATTACH:
      g_hInst = hInstance;
      dll_process_attach( hInstance, true );
      LoadR();
      break;

    case DLL_PROCESS_DETACH:
      dll_process_detach( hInstance, true );
      break;

    case DLL_THREAD_ATTACH:
      dll_thread_attach( true, true );
      break;

    case DLL_THREAD_DETACH:
      dll_thread_detach( true, true );
      break;
    default:
      break;  
    }
    return true;
  }
}
