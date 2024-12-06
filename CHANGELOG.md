# Changelog

## [0.2.2](https://github.com/fxwiegand/predictosaurus/compare/v0.2.1...v0.2.2) (2024-12-06)


### Bug Fixes

* Improve error messages ([c00a01f](https://github.com/fxwiegand/predictosaurus/commit/c00a01f64d205d53c8837b8381da758c8552199b))
* Improve error messages ([3f41035](https://github.com/fxwiegand/predictosaurus/commit/3f41035a54c39d093eea1c6ee52f2362ec8006df))

## [0.2.1](https://github.com/fxwiegand/predictosaurus/compare/v0.2.0...v0.2.1) (2024-12-06)


### Bug Fixes

* Fix out dir creation ([19dad8d](https://github.com/fxwiegand/predictosaurus/commit/19dad8dfec95fc7d742f645551a4b4c9860b60a3))


### Miscellaneous Chores

* release 0.2.1 ([054d16a](https://github.com/fxwiegand/predictosaurus/commit/054d16af036fb56ec93bb016ca99cc8cc84126fb))

## [0.2.0](https://github.com/fxwiegand/predictosaurus/compare/v0.1.1...v0.2.0) (2024-12-05)


### Features

* Add `--min-prop-present` parameter to build subcommand ([c4542cb](https://github.com/fxwiegand/predictosaurus/commit/c4542cb108d81db203ef1fbb0ab2f55b452ec97f))


### Bug Fixes

* Ignore samples present in calls file that are not present in given observations ([2bff91a](https://github.com/fxwiegand/predictosaurus/commit/2bff91a88314031c99d811b46e87d429da52702f))

## [0.1.1](https://github.com/fxwiegand/predictosaurus/compare/v0.1.0...v0.1.1) (2024-12-04)


### Miscellaneous Chores

* release 0.1.1 ([8551a11](https://github.com/fxwiegand/predictosaurus/commit/8551a110d543be8b84acdfe0f0e0e2692f17ec15))

## 0.1.0 (2024-12-04)


### Features

* Add reverse paths for handling features on reverse strand ([897bd81](https://github.com/fxwiegand/predictosaurus/commit/897bd81193fd062b7adae954e373660ecb511607))
* Add support for features on backward strand ([e6b2141](https://github.com/fxwiegand/predictosaurus/commit/e6b21411d3fdce1b54a36f115e8cb76d099ea030))
* Allow conversion of haplotype path to protein ([3fa2c75](https://github.com/fxwiegand/predictosaurus/commit/3fa2c7598268cf7f197aa8943f25b992fb80780e))
* Calculate 0-based position of variant on reverse strand ([8321619](https://github.com/fxwiegand/predictosaurus/commit/8321619843be637c57aad43d671701a93d85b971))
* Enable serialization and deserialization for Graph ([200b4dd](https://github.com/fxwiegand/predictosaurus/commit/200b4dd8cba77efcd3aef66c0151d0bc845c7488))
* Implement show command ([ba4ea1d](https://github.com/fxwiegand/predictosaurus/commit/ba4ea1d1776fdab0d180a9cb96d5b1c62fca9122))
* Implement show command ([a2f243e](https://github.com/fxwiegand/predictosaurus/commit/a2f243ecca5cafd9f65724e103acece84b4bfac8))
* Implement show command ([3bc8d30](https://github.com/fxwiegand/predictosaurus/commit/3bc8d3080dae4e5106a9b7113e08fa92c514def3))
* Implement show command for HTML format ([4af776b](https://github.com/fxwiegand/predictosaurus/commit/4af776b69a00b41a73ec8d8e874195d61be291cd))
* Implement show command for tsv format ([87be0d9](https://github.com/fxwiegand/predictosaurus/commit/87be0d984711520da761e534b540fca3f7a3b8f6))
* Implement show command for tsv format ([4616bc6](https://github.com/fxwiegand/predictosaurus/commit/4616bc6b70c4df1916f8ff68301fe5f2756ee249))
* Implement show command for tsv format ([f0a8965](https://github.com/fxwiegand/predictosaurus/commit/f0a89655851eaf15628d247d6a5857c825e4f87a))
* Restructure into multiple subcommands ([#13](https://github.com/fxwiegand/predictosaurus/issues/13)) ([4462a23](https://github.com/fxwiegand/predictosaurus/commit/4462a233f037cdccc87ec5acfc40af5f0f69d458))
* Support insertions longer 3 ([252d1d6](https://github.com/fxwiegand/predictosaurus/commit/252d1d64577ab4046b0edc4a2795ba3e2d970932))
* Test graph serialization using serde_json ([5074b49](https://github.com/fxwiegand/predictosaurus/commit/5074b49001a8228c67a3c645ae73469812823cec))
* Use duckdb to write graphs ([#35](https://github.com/fxwiegand/predictosaurus/issues/35)) ([19955ea](https://github.com/fxwiegand/predictosaurus/commit/19955eadaac9d3c0f9dadb7480e8265dca71a285))


### Bug Fixes

* Create path table in bd before writing to it ([c335418](https://github.com/fxwiegand/predictosaurus/commit/c3354189f046ce9a24cd78baa5b5349435a790fe))
* Fix impact color scheme ([b8bed5e](https://github.com/fxwiegand/predictosaurus/commit/b8bed5e733d49863efead1100ad1bb18c15755b9))
* Fix impact for amino acid altering variants ([38548a0](https://github.com/fxwiegand/predictosaurus/commit/38548a0e35ecfa95db48069a0c36012b61365abb))
* Handle variants at start and end of reference ([794a7e0](https://github.com/fxwiegand/predictosaurus/commit/794a7e081eb1775cb1b9c83639c080980efda6b3))
* Refactor field names ([2378127](https://github.com/fxwiegand/predictosaurus/commit/2378127b45be953b96f2f61aef799bf6a8e2f6fa))
* Refine behaviour for invalid EventProbs ([c1a300d](https://github.com/fxwiegand/predictosaurus/commit/c1a300d7ba7f4462e2774a869e31befd0d945c3b))
* Refine behaviour for invalid EventProbs ([4342223](https://github.com/fxwiegand/predictosaurus/commit/4342223da59defd6eb39c76fce4225765dc3efa9))
* Refine behaviour for invalid EventProbs ([ac541b1](https://github.com/fxwiegand/predictosaurus/commit/ac541b1382f5304913c8aae762e663a96e63151d))
* Refine behaviour for invalid EventProbs ([d06b7c4](https://github.com/fxwiegand/predictosaurus/commit/d06b7c45849c604f7dbb9704ddc89290b9c08620))
* Refine behaviour for invalid EventProbs ([89ae172](https://github.com/fxwiegand/predictosaurus/commit/89ae172f7948ff94469d9a0f7dac4b6615690f72))
* Remove duplicate ref nodes for multiple variants with same position ([daf9ee0](https://github.com/fxwiegand/predictosaurus/commit/daf9ee0718050cd4cf9782471f82fa036fd05408))
* Use json5 to serialize and deserialize EventProbs ([3330305](https://github.com/fxwiegand/predictosaurus/commit/333030533326cb24acaec31ff776014a60b7cbc9))
