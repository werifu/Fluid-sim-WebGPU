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

export function showFPS(fps: number) {
  let fpsDiv = document.querySelector<HTMLDivElement>('#fps');
  if (!fpsDiv) {
    fpsDiv = document.createElement('div');
    document.body.appendChild(fpsDiv);
  }
  fpsDiv.id = 'fps';
  fpsDiv.style.position = 'fixed';
  fpsDiv.style.top = '0';
  fpsDiv.style.left = '0';

  fpsDiv.innerHTML = 'FPS: ' + fps.toFixed(1);
  if (fps > 45) {
    fpsDiv.style.color = 'green';
  } else if (fps >= 30 && fps <= 45) {
    fpsDiv.style.color = 'yellow';
  } else {
    fpsDiv.style.color = 'red';
  }
}