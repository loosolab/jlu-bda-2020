import { ComponentFixture, TestBed } from '@angular/core/testing';

import { GraphTfComponent } from './graph-tf.component';

describe('GraphTfComponent', () => {
  let component: GraphTfComponent;
  let fixture: ComponentFixture<GraphTfComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ GraphTfComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(GraphTfComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
