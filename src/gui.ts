import { GUI } from 'dat.gui';
import { Scene } from './scene';
import { FlipFluid } from './flipFluid';

export function setupGUI(scene: Scene) {
  const gui = new GUI({ width: 300 });
  gui.add(scene, 'compensateDrift');
  gui.add(scene, 'separateParticles');
  // when we change the number of particles, we need to refresh the fluid
  gui.add(scene, 'dynNumX', 0, 200).onFinishChange((_) => {
    const {
      density,
      simWidth,
      simHeight,
      spacing,
      particleRadius,
      dynNumX,
      dynNumY,
    } = scene;
    scene.fluid = new FlipFluid({ density, simWidth, simHeight, spacing, particleRadius, dynNumX, dynNumY });
  });
  gui.add(scene, 'dynNumY', 0, 200).onFinishChange((_) => {
    const {
      density,
      simWidth,
      simHeight,
      spacing,
      particleRadius,
      dynNumX,
      dynNumY,
    } = scene;
    scene.fluid = new FlipFluid({ density, simWidth, simHeight, spacing, particleRadius, dynNumX, dynNumY });
  });
  // gui.add(scene, 'dynNumY', 0, 200);
}