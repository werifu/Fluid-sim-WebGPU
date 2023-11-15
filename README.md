# Fluid-sim-WebGPU
Final project of Advanced Digital Design

## Development Settings

1. Install LiveServer extension in vscode.
2. Add the following settings to `.vscode/settings.json` in the workspace directory.

*DO NOT USE THIS CERTIFICATE IN PRODUCTION ENVIRONMENT*

```json
{
  "liveServer.settings.https": {
    "enable": true, //set it true to enable the feature.
    "cert": "Your-path/Fluid-sim-WebGPU/dev/RootCA.crt", //full path of the certificate
    "key": "Your-path/Fluid-sim-WebGPU/dev/RootCA.key", //full path of the private key
    "passphrase": "12345"
  }
}
```
