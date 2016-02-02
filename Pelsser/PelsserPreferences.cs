using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DVPLI;
using Mono.Addins;
namespace Pelsser
{

    [SettingsContainer("Pelsser Squared Gaussian Model", SettingType = SettingType.Valuation)]
    [Extension("/Fairmat/UserSettings")]
    [Serializable]
    public class PelsserPreferences: ISettings
    {
        [SettingDescription("Negative Rates Correction")]
        public bool NegativeRatesCorrection = false;
    }
}
